"""
Load VRS-annotated VCF files into pipeline tables.

Usage:
    python manage.py load_vcf [--vcf-out /path/to/vcf_out]

Sources (in priority order for Variant fields):
  1. ENIGMA          → variant, variant_genomic_coordinates, variant_enigma
  2. ClinVar         → variant (upsert), variant_genomic_coordinates, variant_clinvar, report_clinvar
  3. LOVD            → variant (upsert), variant_genomic_coordinates, variant_lovd, report_lovd
  4. exLOVD          → variant (upsert), variant_genomic_coordinates, variant_exlovd
  5. gnomAD v2       → variant, variant_genomic_coordinates, variant_gnomad, report_gnomad (data_type=exome)
  6. gnomAD v3       → variant, variant_genomic_coordinates, variant_gnomad, report_gnomad (data_type=genome)
  7. gnomAD v4       → variant, variant_genomic_coordinates, variant_gnomad, report_gnomad (data_type=joint)
  8. gnomAD v4.1 joint  → variant, variant_genomic_coordinates, variant_gnomad, report_gnomad (data_type=joint)
  9. gnomAD v4.1 exome  → variant, variant_genomic_coordinates, variant_gnomad, report_gnomad (data_type=exome)
"""

import csv
import pickle
import re
from pathlib import Path

import pysam
from django.core.management.base import BaseCommand
from django.db import connections, transaction

from data.models import (
    Genomic_Coordinates,
    Report_in_ClinVar,
    Report_in_GnomAD,
    Report_in_LOVD,
    Variant,
    Variant_in_ClinVar,
    Variant_in_ENIGMA,
    Variant_in_ExLOVD,
    Variant_in_GnomAD,
    Variant_in_LOVD,
    Variant_in_Other,
)

DB = 'pipeline'


# ---------------------------------------------------------------------------
# Helpers (identical logic to load_vcf_to_db.py)
# ---------------------------------------------------------------------------

def info_str(rec, key, default='-'):
    """
    Any INFO field declared Number='.' in these VCF headers (which is all of
    them) gets split by pysam on every literal comma in its value, e.g.
    'lc_lof, os_lof' -> ('lc_lof', ' os_lof') or free text like
    'Burke (Brisbane,AU)' -> ('...Brisbane', 'AU)...'). Rejoining with a bare
    ',' (no added/stripped whitespace) exactly reverses str.split(','),
    whether or not the original had a space after the comma, instead of
    corrupting the text.
    """
    try:
        val = rec.info.get(key)
    except (KeyError, ValueError):
        return default
    if val is None:
        return default
    if isinstance(val, (tuple, list)):
        # A bare INFO key with no '=value' (e.g. 'frequency' alone in the
        # record) parses as an empty tuple rather than None; treat that the
        # same as absent.
        return ','.join(str(v) for v in val).rstrip('\r') if val else default
    return str(val).rstrip('\r')


def alt_digest(rec):
    ids = rec.info.get('VRS_Allele_IDs')
    if ids and len(ids) >= 2:
        return ids[1]
    return None


def alt_vrs_id(rec):
    ids = rec.info.get('VRS_Allele_IDs')
    if ids and len(ids) >= 2:
        return ids[1]
    return None


def populations_by_stat_infix(rec, infix, exclude=frozenset()):
    """
    Build {population: {ac, an, af, ac_hom}} from INFO fields named
    AC_<infix>_<population> / AN_<infix>_<population> / AF_<infix>_<population> /
    nhomalt_<infix>_<population> (infix may be '' when the VCF has no data-type
    infix, e.g. an exome-only release). Used for gnomAD v4/v4.1 joint and exome.
    """
    prefix = f'AC_{infix}_' if infix else 'AC_'
    pattern = re.compile(rf'^{re.escape(prefix)}([a-z]+)$')
    pops = {}
    for key in rec.info.keys():
        m = pattern.match(key)
        if not m or m.group(1) in exclude:
            continue
        pop = m.group(1)
        an_key  = f'AN_{infix}_{pop}'      if infix else f'AN_{pop}'
        af_key  = f'AF_{infix}_{pop}'      if infix else f'AF_{pop}'
        hom_key = f'nhomalt_{infix}_{pop}' if infix else f'nhomalt_{pop}'
        pops[pop] = {
            'ac':     info_str(rec, key),
            'an':     info_str(rec, an_key),
            'af':     info_str(rec, af_key),
            'ac_hom': info_str(rec, hom_key),
        }
    return pops


def populations_by_dataset(rec, dataset, exclude=frozenset()):
    """
    Build {population: {ac, an, af, ac_hom}} from INFO fields named
    <dataset>_<POPULATION>_ac / _an / _af / _ac_hom. Used for gnomAD v2/v3,
    where the population code prefixes the stat rather than the reverse.
    """
    pattern = re.compile(rf'^{dataset}_([A-Z]+)_ac$')
    pops = {}
    for key in rec.info.keys():
        m = pattern.match(key)
        if not m or m.group(1) in exclude:
            continue
        pop_upper = m.group(1)
        pops[pop_upper.lower()] = {
            'ac':     info_str(rec, f'{dataset}_{pop_upper}_ac'),
            'an':     info_str(rec, f'{dataset}_{pop_upper}_an'),
            'af':     info_str(rec, f'{dataset}_{pop_upper}_af'),
            'ac_hom': info_str(rec, f'{dataset}_{pop_upper}_ac_hom'),
        }
    return pops


_CHR_TO_NC = {
    '13': 'NC_000013.11', 'chr13': 'NC_000013.11',
    '17': 'NC_000017.11', 'chr17': 'NC_000017.11',
}


def bare_chrom(chrom):
    return chrom[3:] if chrom.startswith('chr') else chrom


def genomic_hgvs(rec):
    """Construct NC_XXXXXX.XX:g.XXXX genomic HGVS from a pysam VCF record."""
    nc = _CHR_TO_NC.get(rec.chrom)
    if not nc:
        return '-'
    pos = rec.pos + 1  # pysam is 0-based; HGVS uses 1-based positions
    ref = rec.ref
    alt = rec.alts[0] if rec.alts else None
    if alt is None:
        return '-'
    if len(ref) == 1 and len(alt) == 1:
        return f'{nc}:g.{pos}{ref}>{alt}'
    if len(ref) > 1 and len(alt) == 1:
        # deletion with VCF anchor base at pos; deleted bases start at pos+1
        if len(ref) == 2:
            return f'{nc}:g.{pos + 1}del'
        return f'{nc}:g.{pos + 1}_{pos + len(ref) - 1}del'
    if len(ref) == 1 and len(alt) > 1:
        # insertion with VCF anchor base at pos
        return f'{nc}:g.{pos}_{pos + 1}ins{alt[1:]}'
    # complex indel — represent as substitution
    return f'{nc}:g.{pos}_{pos + len(ref) - 1}delins{alt}'


def load_pickle(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


# Stars per ClinVar review status, for choosing between VCVs that
# resolve to the same VRS digest
CLINVAR_REVIEW_STATUS_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, conflicting classifications': 1,
    'criteria provided, single submitter': 1,
}


def clinvar_aggregate_fields(rec):
    """variant_clinvar fields from the VCV_ INFO columns written by clinVarParse.py"""
    variation_id = info_str(rec, 'VCV_VariationID')
    return {
        'Source_URL': (f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/'
                       if variation_id != '-' else '-'),
        'Variation_ID':           variation_id,
        'VCV_Accession':          info_str(rec, 'VCV_Accession'),
        'VCV_Version':            info_str(rec, 'VCV_Version'),
        'Clinical_Significance':  info_str(rec, 'VCV_ClinicalSignificance'),
        'Review_Status':          info_str(rec, 'VCV_ReviewStatus'),
        'Date_Last_Evaluated':    info_str(rec, 'VCV_DateLastEvaluated'),
        'Number_Of_Submissions':  info_str(rec, 'VCV_NumberOfSubmissions'),
        'Number_Of_Submitters':   info_str(rec, 'VCV_NumberOfSubmitters'),
        'Most_Recent_Submission': info_str(rec, 'VCV_MostRecentSubmission'),
        'Explanation':            info_str(rec, 'VCV_Explanation'),
        'Conditions':             info_str(rec, 'VCV_Conditions'),
        'Condition_DB_IDs':       info_str(rec, 'VCV_ConditionDB_IDs'),
        'Date_Last_Updated':      info_str(rec, 'VCV_DateLastUpdated'),
    }


def clinvar_aggregate_rank(fields):
    """Higher is better: review stars, then most recent evaluation"""
    date = fields['Date_Last_Evaluated']
    return (CLINVAR_REVIEW_STATUS_STARS.get(fields['Review_Status'], 0),
            date if date != '-' else '')


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------

class Command(BaseCommand):
    help = 'Load VRS-annotated VCF files into pipeline tables'

    def add_arguments(self, parser):
        parser.add_argument('--gene-config', required=True,
                            help='Path to gene config CSV (must have a "symbol" column)')
        for source in ('enigma', 'clinvar', 'lovd', 'exlovd',
                       'gnomad-v2', 'gnomad-v3',
                       'gnomad-v41-joint', 'gnomad-v41-exome',
                       'functional-assay'):
            parser.add_argument(f'--{source}-vcf',  required=True, help=f'Path to {source} VCF (.vcf.gz)')
            parser.add_argument(f'--{source}-pkl',  required=True, help=f'Path to {source} pickle (.dicts.pkl)')

    def handle(self, *args, **options):
        with open(options['gene_config'], newline='') as f:
            reader = csv.DictReader(f)
            self._known_symbols = {row['symbol'] for row in reader}

        def paths(source):
            key = source.replace('-', '_')
            return Path(options[f'{key}_vcf']), Path(options[f'{key}_pkl'])

        loaders = [
            ('ENIGMA',               self._load_enigma,              *paths('enigma')),
            ('ClinVar',              self._load_clinvar,             *paths('clinvar')),
            ('LOVD',                 self._load_lovd,                *paths('lovd')),
            ('exLOVD',               self._load_exlovd,              *paths('exlovd')),
            ('gnomAD v2',            self._load_gnomad_v2,           *paths('gnomad-v2')),
            ('gnomAD v3',            self._load_gnomad_v3,           *paths('gnomad-v3')),
            ('gnomAD v4.1 joint',    self._load_gnomad_v41_joint,    *paths('gnomad-v41-joint')),
            ('gnomAD v4.1 exome',    self._load_gnomad_v41_exome,    *paths('gnomad-v41-exome')),
            ('functional assays',    self._load_functional_assays,   *paths('functional-assay')),
        ]

        self._flush()

        for label, fn, vcf_path, pkl_path in loaders:
            self.stdout.write(f'Loading {label} ({vcf_path.name}) ...', ending=' ')
            pkl = load_pickle(pkl_path)
            with transaction.atomic(using=DB):
                n = fn(vcf_path, pkl)
            self.stdout.write(str(n))

        self._print_summary()

    def _filter_gene_symbol(self, raw):
        """Return the first symbol in a comma-separated list that is in the gene config.
        Falls back to raw if none match (e.g. a source we don't recognize)."""
        if not raw or raw == '-':
            return raw
        for part in raw.split(','):
            part = part.strip()
            if part in self._known_symbols:
                return part
        return raw

    # ------------------------------------------------------------------
    # Flush
    # ------------------------------------------------------------------

    def _flush(self):
        """Clear all pipeline tables via TRUNCATE CASCADE before reloading."""
        self.stdout.write('Flushing existing data ...')
        tables = [
            'report_gnomad', 'variant_gnomad',
            'variant_exlovd',
            'report_lovd', 'variant_lovd',
            'report_clinvar', 'variant_clinvar',
            'variant_enigma',
            'variant_other',
            'variant_genomic_coordinates',
            'variant',
        ]
        with connections[DB].cursor() as cur:
            for table in tables:
                cur.execute(f'TRUNCATE TABLE "{table}" CASCADE')

    # ------------------------------------------------------------------
    # Shared ORM helpers
    # ------------------------------------------------------------------

    def _upsert_variant(self, d, vrs_data, **fields):
        """Create variant; on conflict fill only fields still at default '-'."""
        variant, created = Variant.objects.using(DB).get_or_create(
            VRS_Digest=d,
            defaults={**fields, 'VRS': vrs_data, 'CA_ID': None},
        )
        if not created:
            update_fields = []
            for attr, val in fields.items():
                if getattr(variant, attr) == '-' and val != '-':
                    setattr(variant, attr, val)
                    update_fields.append(attr)
            if vrs_data and variant.VRS is None:
                variant.VRS = vrs_data
                update_fields.append('VRS')
            if update_fields:
                variant.save(using=DB, update_fields=update_fields)
        return variant

    def _upsert_coords(self, variant, **fields):
        Genomic_Coordinates.objects.using(DB).get_or_create(
            VRS_Digest=variant, assembly='GRCh38', defaults=fields,
        )

    # ------------------------------------------------------------------
    # Source loaders
    # ------------------------------------------------------------------

    def _load_enigma(self, vcf_path, pkl):
        vcf = pysam.VariantFile(str(vcf_path))
        n = 0
        for rec in vcf:
            d = alt_digest(rec)
            if not d:
                continue
            vrs_data = pkl.get(alt_vrs_id(rec))

            # ENIGMA is authoritative — always overwrite existing variant fields
            variant, created = Variant.objects.using(DB).get_or_create(
                VRS_Digest=d,
                defaults={
                    'Gene_Symbol':        self._filter_gene_symbol(info_str(rec, 'Gene_symbol')),
                    'Reference_Sequence': info_str(rec, 'Reference_sequence'),
                    'HGVS_cDNA':          info_str(rec, 'HGVS_cDNA'),
                    'BIC_Nomenclature':   info_str(rec, 'BIC_Nomenclature'),
                    'HGVS_Protein':       info_str(rec, 'HGVS_protein'),
                    'Protein_Change':     info_str(rec, 'Abbrev_AA_change'),
                    'CA_ID':              None,
                    'VRS':                vrs_data,
                    'Synonyms':           '-',
                },
            )
            if not created:
                variant.Gene_Symbol        = self._filter_gene_symbol(info_str(rec, 'Gene_symbol'))
                variant.Reference_Sequence = info_str(rec, 'Reference_sequence')
                variant.HGVS_cDNA          = info_str(rec, 'HGVS_cDNA')
                variant.BIC_Nomenclature   = info_str(rec, 'BIC_Nomenclature')
                variant.HGVS_Protein       = info_str(rec, 'HGVS_protein')
                variant.Protein_Change     = info_str(rec, 'Abbrev_AA_change')
                if vrs_data:
                    variant.VRS = vrs_data
                variant.save(using=DB)

            self._upsert_coords(variant,
                hgvs=genomic_hgvs(rec),
                end_pos=str(rec.pos + len(rec.ref) - 1),
                chr=bare_chrom(rec.chrom), pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            Variant_in_ENIGMA.objects.using(DB).create(
                VRS_Digest=variant,
                Condition_ID_type=info_str(rec, 'Condition_ID_type'),
                Condition_ID_value=info_str(rec, 'Condition_ID_value'),
                Condition_category=info_str(rec, 'Condition_category'),
                Clinical_significance=info_str(rec, 'Clinical_significance'),
                Date_last_evaluated=info_str(rec, 'Date_last_evaluated'),
                Assertion_method=info_str(rec, 'Assertion_method'),
                Assertion_method_citation=info_str(rec, 'Assertion_method_citation'),
                Clinical_significance_citations=info_str(rec, 'Clinical_significance_citations'),
                Comment_on_clinical_significance=info_str(rec, 'Comment_on_clinical_significance'),
                Collection_method=info_str(rec, 'Collection_method'),
                Allele_origin=info_str(rec, 'Allele_origin'),
                ClinVarAccession=info_str(rec, 'ClinVarAccession'),
                Pathogenicity=info_str(rec, 'Clinical_significance'),
            )
            n += 1
        vcf.close()
        return n

    def _load_clinvar(self, vcf_path, pkl):
        vcf = pysam.VariantFile(str(vcf_path))
        n = 0
        for rec in vcf:
            d = alt_digest(rec)
            if not d:
                continue
            vrs_data = pkl.get(alt_vrs_id(rec))
            aggregate = clinvar_aggregate_fields(rec)

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol=self._filter_gene_symbol(info_str(rec, 'Symbol')),
                Reference_Sequence='-', HGVS_cDNA='-',
                BIC_Nomenclature=info_str(rec, 'BIC_Nomenclature'),
                HGVS_Protein='-', Protein_Change='-',
                Synonyms=info_str(rec, 'Synonyms'),
            )

            self._upsert_coords(variant,
                hgvs=genomic_hgvs(rec),
                end_pos=str(rec.pos + len(rec.ref) - 1),
                chr=bare_chrom(rec.chrom), pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            clinvar_rec, created = Variant_in_ClinVar.objects.using(DB).get_or_create(
                VRS_Digest=variant, defaults=aggregate,
            )
            if (not created and aggregate['VCV_Accession'] != clinvar_rec.VCV_Accession
                    and clinvar_aggregate_rank(aggregate) > clinvar_aggregate_rank(
                        {k: getattr(clinvar_rec, k) for k in aggregate})):
                for k, v in aggregate.items():
                    setattr(clinvar_rec, k, v)
                clinvar_rec.save(using=DB)

            Report_in_ClinVar.objects.using(DB).create(
                VRS_Digest=clinvar_rec,
                Clinical_Significance=info_str(rec, 'ClinicalSignificance'),
                Date_Last_Updated=info_str(rec, 'DateLastUpdated'),
                DateSignificanceLastEvaluated=info_str(rec, 'DateSignificanceLastEvaluated'),
                Submitter=info_str(rec, 'Submitter'),
                SCV=info_str(rec, 'SCV'),
                SCV_Version=info_str(rec, 'SCV_Version'),
                Allele_Origin=info_str(rec, 'Origin'),
                Method=info_str(rec, 'Method'),
                Description=info_str(rec, 'Description'),
                Summary_Evidence=info_str(rec, 'SummaryEvidence'),
                Review_Status=info_str(rec, 'ReviewStatus'),
                Condition_Type=info_str(rec, 'ConditionType'),
                Condition_Value=info_str(rec, 'ConditionValue'),
                Condition_DB_ID=info_str(rec, 'ConditionDB_ID'),
            )
            n += 1
        vcf.close()
        return n

    def _load_lovd(self, vcf_path, pkl):
        vcf = pysam.VariantFile(str(vcf_path))
        n = 0
        for rec in vcf:
            d = alt_digest(rec)
            if not d:
                continue
            vrs_data = pkl.get(alt_vrs_id(rec))
            dbid = info_str(rec, 'DBID')
            source_url = (
                f'https://databases.lovd.nl/shared/variants/{dbid}'
                if dbid != '-' else '-'
            )

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol=self._filter_gene_symbol(info_str(rec, 'geneid')),
                Reference_Sequence='-',
                HGVS_cDNA=info_str(rec, 'cDNA'),
                BIC_Nomenclature='-', HGVS_Protein='-', Protein_Change='-',
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=genomic_hgvs(rec),
                end_pos=str(rec.pos + len(rec.ref) - 1),
                chr=bare_chrom(rec.chrom), pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            lovd_rec, _ = Variant_in_LOVD.objects.using(DB).get_or_create(
                VRS_Digest=variant,
                defaults={
                    'Source_URL': source_url,
                    'Variant_haplotype': info_str(rec, 'haplotype'),
                },
            )

            Report_in_LOVD.objects.using(DB).create(
                VRS_Digest=lovd_rec,
                Variant_frequency=info_str(rec, 'frequency'),
                Individuals=info_str(rec, 'individuals'),
                Variant_effect=info_str(rec, 'variant_effect'),
                Genetic_origin=info_str(rec, 'genetic_origin'),
                Submitters=info_str(rec, 'submitters'),
                Functional_analysis_technique=info_str(rec, 'functionalanalysis_technique'),
                Functional_analysis_result=info_str(rec, 'functionalanalysis_result'),
                Created_date=info_str(rec, 'created_date'),
                Edited_date=info_str(rec, 'edited_date'),
                DBID=dbid,
                Remarks=info_str(rec, 'remarks'),
                Classification=info_str(rec, 'classification'),
                Submission_ID=info_str(rec, 'submission_id'),
            )
            n += 1
        vcf.close()
        return n

    def _load_exlovd(self, vcf_path, pkl):
        vcf = pysam.VariantFile(str(vcf_path))
        n = 0
        for rec in vcf:
            d = alt_digest(rec)
            if not d:
                continue
            vrs_data = pkl.get(alt_vrs_id(rec))

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol='-', Reference_Sequence='-',
                HGVS_cDNA=info_str(rec, 'dna_change'),
                BIC_Nomenclature=info_str(rec, 'bic_dna_change'),
                HGVS_Protein='-',
                Protein_Change=info_str(rec, 'protein_change'),
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=genomic_hgvs(rec),
                end_pos=str(rec.pos + len(rec.ref) - 1),
                chr=bare_chrom(rec.chrom), pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            Variant_in_ExLOVD.objects.using(DB).get_or_create(
                VRS_Digest=variant,
                defaults={
                    'Posterior_P':               info_str(rec, 'posterior_p'),
                    'IARC_Class':                info_str(rec, 'iarc_class'),
                    'Missense_Analysis_Prior_P': info_str(rec, 'missense_analysis_prior_p'),
                    'Combined_Prior_P':          info_str(rec, 'combined_prior_p'),
                    'Segregation_LR':            info_str(rec, 'segregation_lr'),
                    'Pathology_LR':              info_str(rec, 'pathology_lr') or None,
                    'Co_Occurrence_LR':          info_str(rec, 'co_occurrence_lr'),
                    'Case_Control_LR':           info_str(rec, 'case_control_lr') or None,
                    'Comments':                  info_str(rec, 'key_observational_reference'),
                },
            )
            n += 1
        vcf.close()
        return n

    def _load_gnomad(self, vcf_path, pkl, version, data_type, ac_key, an_key, af_key,
                     variant_id_key, faf_key, faf_pop_key, build_populations,
                     coverage_key=None):
        vcf = pysam.VariantFile(str(vcf_path))
        n = 0
        for rec in vcf:
            d = alt_digest(rec)
            if not d:
                continue
            vrs_data = pkl.get(alt_vrs_id(rec))

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol='-', Reference_Sequence='-', HGVS_cDNA='-',
                BIC_Nomenclature='-', HGVS_Protein='-', Protein_Change='-',
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=genomic_hgvs(rec),
                end_pos=str(rec.pos + len(rec.ref) - 1),
                chr=bare_chrom(rec.chrom), pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            gnomad_rec, _ = Variant_in_GnomAD.objects.using(DB).get_or_create(
                VRS_Digest=variant, defaults={'Source_URL': '-'},
            )

            Report_in_GnomAD.objects.using(DB).get_or_create(
                VRS_Digest=gnomad_rec,
                version=version,
                data_type=data_type,
                defaults={
                    'Variant_id':              info_str(rec, variant_id_key),
                    'Flags':                   info_str(rec, 'flags'),
                    'coverage':                info_str(rec, coverage_key) if coverage_key else '-',
                    'Allele_count':            info_str(rec, ac_key),
                    'Allele_number':           info_str(rec, an_key),
                    'Allele_frequency':        info_str(rec, af_key),
                    'faf95_popmax':            info_str(rec, faf_key),
                    'faf95_popmax_population': info_str(rec, faf_pop_key),
                    'populations':             build_populations(rec),
                },
            )
            n += 1
        vcf.close()
        return n

    def _load_gnomad_v2(self, vcf_path, pkl):
        def populations(rec):
            return populations_by_dataset(rec, 'genome', exclude={'FEMALE', 'MALE'})
        return self._load_gnomad(vcf_path, pkl, version='v2', data_type='exome',
            ac_key='genome_ac', an_key='genome_an', af_key='genome_af',
            variant_id_key='variantId',
            faf_key='genome_popmax', faf_pop_key='genome_popmax_population',
            build_populations=populations,
        )

    def _load_gnomad_v3(self, vcf_path, pkl):
        def populations(rec):
            return populations_by_dataset(rec, 'genome', exclude={'FEMALE', 'MALE'})
        return self._load_gnomad(vcf_path, pkl, version='v3', data_type='genome',
            ac_key='genome_ac', an_key='genome_an', af_key='genome_af',
            variant_id_key='variant_id',
            faf_key='genome_popmax', faf_pop_key='genome_popmax_population',
            build_populations=populations, coverage_key='coverage',
        )

    def _load_gnomad_v4(self, vcf_path, pkl):
        def populations(rec):
            return populations_by_stat_infix(rec, 'joint', exclude={'raw'})
        return self._load_gnomad(vcf_path, pkl, version='v4', data_type='joint',
            ac_key='AC_joint', an_key='AN_joint', af_key='AF_joint',
            variant_id_key='variant_id',
            faf_key='fafmax_faf95_max_joint', faf_pop_key='fafmax_faf95_max_gen_anc_joint',
            build_populations=populations,
        )

    def _load_gnomad_v41_joint(self, vcf_path, pkl):
        def populations(rec):
            return populations_by_stat_infix(rec, 'joint', exclude={'raw'})
        return self._load_gnomad(vcf_path, pkl, version='v4.1', data_type='joint',
            ac_key='AC_joint', an_key='AN_joint', af_key='AF_joint',
            variant_id_key='variant_id',
            faf_key='fafmax_faf95_max_joint', faf_pop_key='fafmax_faf95_max_gen_anc_joint',
            build_populations=populations, coverage_key='coverage',
        )

    def _load_gnomad_v41_exome(self, vcf_path, pkl):
        def populations(rec):
            return populations_by_stat_infix(rec, '', exclude={'raw', 'grpmax'})
        return self._load_gnomad(vcf_path, pkl, version='v4.1', data_type='exome',
            ac_key='AC', an_key='AN', af_key='AF',
            variant_id_key='variant_id',
            faf_key='fafmax_faf95_max', faf_pop_key='fafmax_faf95_max_gen_anc',
            build_populations=populations, coverage_key='coverage',
        )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def _print_summary(self):
        self.stdout.write('')
        rows = [
            ('variant',              Variant),
            ('variant_genomic_coordinates',  Genomic_Coordinates),
            ('variant_enigma',       Variant_in_ENIGMA),
            ('variant_clinvar',      Variant_in_ClinVar),
            ('report_clinvar',       Report_in_ClinVar),
            ('variant_lovd',         Variant_in_LOVD),
            ('report_lovd',          Report_in_LOVD),
            ('variant_exlovd',       Variant_in_ExLOVD),
            ('variant_gnomad',       Variant_in_GnomAD),
            ('report_gnomad',        Report_in_GnomAD),
            ('variant_other',        Variant_in_Other),
        ]
        for label, model in rows:
            count = model.objects.using(DB).count()
            self.stdout.write(f'{label}: {count}')

    # ------------------------------------------------------------------
    # Functional assay loader
    # ------------------------------------------------------------------

    _VRS_INFO_KEYS = {'VRS_Allele_IDs', 'VRS_Error'}

    def _load_functional_assays(self, vcf_path, pkl):
        vcf = pysam.VariantFile(str(vcf_path))
        n = 0
        for rec in vcf:
            d = alt_digest(rec)
            if not d:
                continue
            vrs_data = pkl.get(alt_vrs_id(rec))

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol=self._filter_gene_symbol(info_str(rec, 'Gene')),
                Reference_Sequence='-',
                HGVS_cDNA=info_str(rec, 'HGVS_Nucleotide_Variant'),
                BIC_Nomenclature='-', HGVS_Protein='-', Protein_Change='-',
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=genomic_hgvs(rec),
                end_pos=str(rec.pos + len(rec.ref) - 1),
                chr=bare_chrom(rec.chrom), pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            data = {}
            for key in rec.info.keys():
                if key in self._VRS_INFO_KEYS:
                    continue
                val = rec.info.get(key)
                if val is None or (isinstance(val, tuple) and len(val) == 0):
                    continue
                if isinstance(val, tuple):
                    data[key] = val[0] if len(val) == 1 else ','.join(str(v) for v in val)
                else:
                    data[key] = str(val)

            Variant_in_Other.objects.using(DB).update_or_create(
                VRS_Digest=variant,
                defaults={'data_type': 'functional_assay_results', 'variant_data': data},
            )
            n += 1
        vcf.close()
        return n
