"""
Load varaico-mined BRCA1/BRCA2 literature mentions into paper / variant_in_paper.

Usage:
    python manage.py add_varaico_papers --raw-tsv <path> --vrs-annotated-vcf <path>

Only links variants that already exist in `variant` (VRS_Digest) — never creates new Variant
rows. Fed by the Luigi chain in workflow/variant_assembly.py:
ExtractVaraicoBRCARegions -> FilterVaraicoToGeneBoundaries -> SortVaraicoVCF ->
VRSAnnotateVaraico -> LoadVaraicoPapersToDatabase (this command).
"""

import csv
import re

import pysam
from django.core.management.base import BaseCommand
from django.db import transaction

from data.models import Paper, Variant, Variant_in_Paper

DB = 'pipeline'

# variantSnippets/variantOrigStrs come from varaico with "<BR/><BR/>"-separated snippets and
# <b>...</b> markup around the matched term.
_SNIPPET_SPLIT_RE = re.compile(r'<BR/>\s*<BR/>')
_HTML_TAG_RE = re.compile(r'<[^>]+>')


def alt_digest(rec):
    ids = rec.info.get('VRS_Allele_IDs')
    if ids and len(ids) >= 2:
        return ids[1]
    return None


def clean_snippets(raw):
    """Split on snippet separators and strip HTML markup, dropping empty pieces."""
    if not raw or raw == '-':
        return []
    pieces = (_HTML_TAG_RE.sub('', p).strip() for p in _SNIPPET_SPLIT_RE.split(raw))
    return [p for p in pieces if p]


class Command(BaseCommand):
    help = 'Load varaico literature-mining data into paper / variant_in_paper'

    def add_arguments(self, parser):
        parser.add_argument('--raw-tsv', required=True,
                            help='Path to varaico_mentions_filtered.tsv (FilterVaraicoToGeneBoundaries output)')
        parser.add_argument('--vrs-annotated-vcf', required=True,
                            help='Path to VRS-annotated varaico VCF (VRSAnnotateVaraico output)')

    def handle(self, *args, **options):
        digest_by_key = self._load_digests(options['vrs_annotated_vcf'])

        n_rows = n_no_digest = n_no_variant = 0
        n_papers_created = n_links_created = n_links_updated = 0

        with open(options['raw_tsv'], newline='') as f, transaction.atomic(using=DB):
            reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader:
                n_rows += 1
                pos = int(row['chromStart']) + 1
                key = f"{row['chrom']}:{pos}:{row['ref']}:{row['alt']}"
                digest = digest_by_key.get(key)
                if not digest:
                    n_no_digest += 1
                    continue

                variant = Variant.objects.using(DB).filter(VRS_Digest=digest).first()
                if variant is None:
                    n_no_variant += 1
                    continue

                paper, created = Paper.objects.using(DB).get_or_create(
                    PMID=row['selectedPmid'],
                    defaults={
                        'Title':   row['title'],
                        'Author':  row['author'],
                        'Year':    row['year'],
                        'Journal': row['journal'],
                        'DOI':     row['doi'],
                    },
                )
                if created:
                    n_papers_created += 1

                mentions = clean_snippets(row['variantSnippets'])
                mentioned_as = clean_snippets(row['variantOrigStrs'])

                vip, created = Variant_in_Paper.objects.using(DB).get_or_create(
                    VRS_Digest=variant, Paper=paper,
                    defaults={'mentions': mentions, 'variant_mentioned_as': mentioned_as},
                )
                if created:
                    n_links_created += 1
                else:
                    new_mentions = [m for m in mentions if m not in vip.mentions]
                    new_mentioned_as = [m for m in mentioned_as if m not in vip.variant_mentioned_as]
                    if new_mentions or new_mentioned_as:
                        vip.mentions += new_mentions
                        vip.variant_mentioned_as += new_mentioned_as
                        vip.save(using=DB, update_fields=['mentions', 'variant_mentioned_as'])
                        n_links_updated += 1

        self.stdout.write(
            f'Rows read: {n_rows}, no VRS digest: {n_no_digest}, '
            f'variant not in DB: {n_no_variant}, papers created: {n_papers_created}, '
            f'links created: {n_links_created}, links updated: {n_links_updated}'
        )

    def _load_digests(self, vcf_path):
        vcf = pysam.VariantFile(vcf_path)
        digest_by_key = {}
        for rec in vcf:
            d = alt_digest(rec)
            if d:
                digest_by_key[rec.id] = d
        vcf.close()
        return digest_by_key
