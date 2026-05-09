"""
Load VRS-annotated VCF files into pipeline tables.

Usage:
    python manage.py load_vcf [--vcf-out /path/to/vcf_out]

Sources (in priority order for Variant fields):
  1. ENIGMA     → variant, genomic_coordinates, variant_enigma
  2. ClinVar    → variant (upsert), genomic_coordinates, variant_clinvar, report_clinvar
  3. LOVD       → variant (upsert), genomic_coordinates, variant_lovd, report_lovd
  4. exLOVD     → variant (upsert), genomic_coordinates, variant_exlovd
  5. gnomAD v2  → variant, genomic_coordinates, variant_gnomad, report_gnomad
  6. gnomAD v3  → variant, genomic_coordinates, variant_gnomad, report_gnomad
  7. gnomAD v4  → variant, genomic_coordinates, variant_gnomad, report_gnomad
"""

import pickle
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
    try:
        val = rec.info.get(key)
    except (KeyError, ValueError):
        return default
    if val is None:
        return default
    if isinstance(val, (tuple, list)):
        return '|'.join(str(v) for v in val)
    return str(val)


def alt_digest(rec):
    ids = rec.info.get('VRS_Allele_IDs')
    if ids and len(ids) >= 2:
        return ids[1].removeprefix('ga4gh:VA.')
    return None


def alt_vrs_id(rec):
    ids = rec.info.get('VRS_Allele_IDs')
    if ids and len(ids) >= 2:
        return ids[1]
    return None


def ucsc_url(chrom, pos):
    c = chrom if chrom.startswith('chr') else f'chr{chrom}'
    return f'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={c}:{pos}-{pos}'


def load_pickle(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------

class Command(BaseCommand):
    help = 'Load VRS-annotated VCF files into pipeline tables'

    def add_arguments(self, parser):
        for source in ('enigma', 'clinvar', 'lovd', 'exlovd',
                       'gnomad-v2', 'gnomad-v3', 'gnomad-v4',
                       'functional-assay'):
            parser.add_argument(f'--{source}-vcf',  required=True, help=f'Path to {source} VCF (.vcf.gz)')
            parser.add_argument(f'--{source}-pkl',  required=True, help=f'Path to {source} pickle (.dicts.pkl)')

    def handle(self, *args, **options):
        def paths(source):
            key = source.replace('-', '_')
            return Path(options[f'{key}_vcf']), Path(options[f'{key}_pkl'])

        loaders = [
            ('ENIGMA',            self._load_enigma,            *paths('enigma')),
            ('ClinVar',           self._load_clinvar,           *paths('clinvar')),
            ('LOVD',              self._load_lovd,              *paths('lovd')),
            ('exLOVD',            self._load_exlovd,            *paths('exlovd')),
            ('gnomAD v2',         self._load_gnomad_v2,         *paths('gnomad-v2')),
            ('gnomAD v3',         self._load_gnomad_v3,         *paths('gnomad-v3')),
            ('gnomAD v4',         self._load_gnomad_v4,         *paths('gnomad-v4')),
            ('functional assays', self._load_functional_assays, *paths('functional-assay')),
        ]

        self._flush()

        for label, fn, vcf_path, pkl_path in loaders:
            self.stdout.write(f'Loading {label} ({vcf_path.name}) ...', ending=' ')
            pkl = load_pickle(pkl_path)
            with transaction.atomic(using=DB):
                n = fn(vcf_path, pkl)
            self.stdout.write(str(n))

        self._print_summary()

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
            'genomic_coordinates',
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
                    'Gene_Symbol':        info_str(rec, 'Gene_symbol'),
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
                variant.Gene_Symbol        = info_str(rec, 'Gene_symbol')
                variant.Reference_Sequence = info_str(rec, 'Reference_sequence')
                variant.HGVS_cDNA          = info_str(rec, 'HGVS_cDNA')
                variant.BIC_Nomenclature   = info_str(rec, 'BIC_Nomenclature')
                variant.HGVS_Protein       = info_str(rec, 'HGVS_protein')
                variant.Protein_Change     = info_str(rec, 'Abbrev_AA_change')
                if vrs_data:
                    variant.VRS = vrs_data
                variant.save(using=DB)

            self._upsert_coords(variant,
                hgvs=info_str(rec, 'HGVS_cDNA'),
                genome_browser_url=ucsc_url(rec.chrom, rec.pos),
                chr=rec.chrom, pos=str(rec.pos),
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
            clinvar_id = info_str(rec, 'ID')
            source_url = (
                f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/'
                if clinvar_id != '-' else '-'
            )

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol=info_str(rec, 'Symbol'),
                Reference_Sequence='-', HGVS_cDNA='-', BIC_Nomenclature='-',
                HGVS_Protein='-', Protein_Change='-',
                Synonyms=info_str(rec, 'Synonyms'),
            )

            self._upsert_coords(variant,
                hgvs=info_str(rec, 'HGVS'),
                genome_browser_url=ucsc_url(rec.chrom, rec.pos),
                chr=rec.chrom, pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            clinvar_rec, _ = Variant_in_ClinVar.objects.using(DB).get_or_create(
                VRS_Digest=variant, defaults={'Source_URL': source_url},
            )

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
                Gene_Symbol=info_str(rec, 'geneid'),
                Reference_Sequence='-',
                HGVS_cDNA=info_str(rec, 'cDNA'),
                BIC_Nomenclature='-', HGVS_Protein='-', Protein_Change='-',
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=info_str(rec, 'gDNA'),
                genome_browser_url=ucsc_url(rec.chrom, rec.pos),
                chr=rec.chrom, pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            lovd_rec, _ = Variant_in_LOVD.objects.using(DB).get_or_create(
                VRS_Digest=variant,
                defaults={
                    'Source_URL': source_url,
                    'Variant_haplotype': info_str(rec, 'haplotype'),
                },
            )

            remarks = rec.info.get('remarks')
            classification = rec.info.get('classification')
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
                Remarks='|'.join(str(v) for v in remarks) if isinstance(remarks, tuple) else (remarks or '-'),
                Classification='|'.join(str(v) for v in classification) if isinstance(classification, tuple) else (classification or '-'),
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
            db_id = info_str(rec, 'db_id')
            source_url = (
                f'https://databases.lovd.nl/shared/variants/{db_id}'
                if db_id != '-' else '-'
            )

            variant = self._upsert_variant(d, vrs_data,
                Gene_Symbol='-', Reference_Sequence='-',
                HGVS_cDNA=info_str(rec, 'dna_change'),
                BIC_Nomenclature=info_str(rec, 'bic_dna_change'),
                HGVS_Protein='-',
                Protein_Change=info_str(rec, 'protein_change'),
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=info_str(rec, 'dna_change'),
                genome_browser_url=ucsc_url(rec.chrom, rec.pos),
                chr=rec.chrom, pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            Variant_in_ExLOVD.objects.using(DB).get_or_create(
                VRS_Digest=variant,
                defaults={
                    'Source_URL':                source_url,
                    'Exon':                      info_str(rec, 'exon'),
                    'DNA_Change':                info_str(rec, 'dna_change'),
                    'BIC_DNA_Change':            info_str(rec, 'bic_dna_change'),
                    'Protein_Change':            info_str(rec, 'protein_change'),
                    'Posterior_P':               info_str(rec, 'posterior_p'),
                    'IARC_Class':                info_str(rec, 'iarc_class'),
                    'DBID':                      db_id,
                    'Missense_Analysis_Prior_P': info_str(rec, 'missense_analysis_prior_p'),
                    'Splicing_Prior_P':          info_str(rec, 'splicing_prior_p'),
                    'Combined_Prior_P':          info_str(rec, 'combined_prior_p'),
                    'Segregation_LR':            info_str(rec, 'segregation_lr'),
                    'Pathology_LR':              info_str(rec, 'pathology_lr') or None,
                    'Co_Occurrence_LR':          info_str(rec, 'co_occurrence_lr'),
                    'Case_Control_LR':           info_str(rec, 'case_control_lr') or None,
                    'Product_Of_LRs':            info_str(rec, 'product_of_lrs'),
                    'Comments':                  info_str(rec, 'comments'),
                },
            )
            n += 1
        vcf.close()
        return n

    def _load_gnomad(self, vcf_path, pkl, version, ac_key, an_key, af_key,
                     variant_id_key, faf_key, faf_pop_key, build_populations):
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
                hgvs=info_str(rec, 'hgvs'),
                genome_browser_url=ucsc_url(rec.chrom, rec.pos),
                chr=rec.chrom, pos=str(rec.pos),
                ref=rec.ref, alt=rec.alts[0] if rec.alts else '-',
            )

            gnomad_rec, _ = Variant_in_GnomAD.objects.using(DB).get_or_create(
                VRS_Digest=variant, defaults={'Source_URL': '-'},
            )

            Report_in_GnomAD.objects.using(DB).get_or_create(
                VRS_Digest=gnomad_rec,
                version=version,
                defaults={
                    'Variant_id':            info_str(rec, variant_id_key),
                    'Flags':                 info_str(rec, 'flags'),
                    'Allele_count':          info_str(rec, ac_key),
                    'Allele_number':         info_str(rec, an_key),
                    'Allele_frequency':      info_str(rec, af_key),
                    'faf95_popmax':          info_str(rec, faf_key),
                    'faf95_popmax_population': info_str(rec, faf_pop_key),
                    'populations':           build_populations(rec),
                },
            )
            n += 1
        vcf.close()
        return n

    def _load_gnomad_v2(self, vcf_path, pkl):
        def populations(rec):
            return {
                'genome_ac': info_str(rec, 'genome_ac'),
                'genome_an': info_str(rec, 'genome_an'),
                'genome_af': info_str(rec, 'genome_af'),
                'exome_ac':  info_str(rec, 'exome_ac'),
                'exome_an':  info_str(rec, 'exome_an'),
                'exome_af':  info_str(rec, 'exome_af'),
            }
        return self._load_gnomad(vcf_path, pkl, version='v2',
            ac_key='genome_ac', an_key='genome_an', af_key='genome_af',
            variant_id_key='variantId',
            faf_key='popmax_genome', faf_pop_key='popmax_population_genome',
            build_populations=populations,
        )

    def _load_gnomad_v3(self, vcf_path, pkl):
        def populations(rec):
            return {
                'genome_ac':     info_str(rec, 'genome_ac'),
                'genome_an':     info_str(rec, 'genome_an'),
                'genome_af':     info_str(rec, 'genome_af'),
                'genome_ac_hom': info_str(rec, 'genome_ac_hom'),
            }
        return self._load_gnomad(vcf_path, pkl, version='v3',
            ac_key='genome_ac', an_key='genome_an', af_key='genome_af',
            variant_id_key='variant_id',
            faf_key='genome_popmax', faf_pop_key='genome_popmax_population',
            build_populations=populations,
        )

    def _load_gnomad_v4(self, vcf_path, pkl):
        def populations(rec):
            return {
                'joint_ac':     info_str(rec, 'AC_joint'),
                'joint_an':     info_str(rec, 'AN_joint'),
                'joint_af':     info_str(rec, 'AF_joint'),
                'genome_ac':    info_str(rec, 'AC_genomes'),
                'genome_an':    info_str(rec, 'AN_genomes'),
                'genome_af':    info_str(rec, 'AF_genomes'),
                'exome_ac':     info_str(rec, 'AC_exomes'),
                'exome_an':     info_str(rec, 'AN_exomes'),
                'exome_af':     info_str(rec, 'AF_exomes'),
                'nhomalt_joint': info_str(rec, 'nhomalt_joint'),
            }
        return self._load_gnomad(vcf_path, pkl, version='v4',
            ac_key='AC_joint', an_key='AN_joint', af_key='AF_joint',
            variant_id_key='variant_id',
            faf_key='faf95_joint', faf_pop_key='grpmax_joint',
            build_populations=populations,
        )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def _print_summary(self):
        self.stdout.write('')
        rows = [
            ('variant',              Variant),
            ('genomic_coordinates',  Genomic_Coordinates),
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
                Gene_Symbol=info_str(rec, 'Gene'),
                Reference_Sequence='-',
                HGVS_cDNA=info_str(rec, 'HGVS_Nucleotide_Variant'),
                BIC_Nomenclature='-', HGVS_Protein='-', Protein_Change='-',
                Synonyms='-',
            )

            self._upsert_coords(variant,
                hgvs=info_str(rec, 'HGVS_Nucleotide_Variant'),
                genome_browser_url=ucsc_url(rec.chrom, rec.pos),
                chr=rec.chrom, pos=str(rec.pos),
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
                    data[key] = val[0] if len(val) == 1 else '|'.join(str(v) for v in val)
                else:
                    data[key] = str(val)

            Variant_in_Other.objects.using(DB).update_or_create(
                VRS_Digest=variant,
                defaults={'data_type': 'functional_assay_results', 'variant_data': data},
            )
            n += 1
        vcf.close()
        return n
