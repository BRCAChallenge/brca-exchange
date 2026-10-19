"""
Fill the VCV-level aggregate classification columns of variant_clinvar, in
place, from a VRS-annotated ClinVar VCF.  Unlike load_vcf it doesn't flush
anything: only existing variant_clinvar rows are updated.

Usage:
    python manage.py load_clinvar_aggregate --clinvar-vcf /path/to/ClinVar.vcf.gz
"""
import pysam
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from data.management.commands.load_vcf import (
    DB, alt_digest, clinvar_aggregate_fields, clinvar_aggregate_rank)
from data.models import Variant_in_ClinVar


class Command(BaseCommand):
    help = 'Fill variant_clinvar aggregate classification columns in place'

    def add_arguments(self, parser):
        parser.add_argument('--clinvar-vcf', required=True,
                            help='Path to the VRS-annotated ClinVar VCF (.vcf.gz)')

    def handle(self, *args, **options):
        best, vcvs_per_digest = {}, {}
        with pysam.VariantFile(options['clinvar_vcf']) as vcf:
            for rec in vcf:
                d = alt_digest(rec)
                if not d:
                    continue
                fields = clinvar_aggregate_fields(rec)
                vcvs_per_digest.setdefault(d, set()).add(fields['VCV_Accession'])
                if d not in best or clinvar_aggregate_rank(fields) > clinvar_aggregate_rank(best[d]):
                    best[d] = fields

        existing = set(Variant_in_ClinVar.objects.using(DB).values_list('VRS_Digest', flat=True))
        missing, extra = set(best) - existing, existing - set(best)
        if missing or extra:
            raise CommandError(f'{len(missing)} VCF digests have no variant_clinvar row and '
                               f'{len(extra)} variant_clinvar rows are not in the VCF; '
                               'nothing written')

        collisions = {d: v for d, v in vcvs_per_digest.items() if len(v) > 1}
        for d, accessions in sorted(collisions.items()):
            self.stdout.write(f'{d}: {sorted(accessions)} -> {best[d]["VCV_Accession"]}')

        field_names = list(next(iter(best.values())).keys())
        rows = [Variant_in_ClinVar(VRS_Digest_id=d, **fields) for d, fields in best.items()]
        with transaction.atomic(using=DB):
            Variant_in_ClinVar.objects.using(DB).bulk_update(rows, field_names, batch_size=1000)
        self.stdout.write(f'Updated {len(rows)} variant_clinvar rows; '
                          f'{len(collisions)} digests had more than one VCV')
