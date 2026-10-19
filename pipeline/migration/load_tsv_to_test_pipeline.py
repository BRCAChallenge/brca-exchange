"""
Load pipeline release TSVs into the test_pipeline schema.

Usage:
    python load_tsv_to_test_pipeline.py \\
        --built-tsv  /path/to/built_with_change_types.tsv \\
        --reports-tsv /path/to/reports_with_change_types.tsv

Pass 1 (built TSV, one row per variant):
    Populates: variant, variant_genomic_coordinates, variant_enigma,
               variant_clinvar, variant_lovd, variant_exlovd, variant_gnomad,
               analysis_priors, analysis_bayesdel,
               analysis_provisional_evidence_codes.

Pass 2 (reports TSV, one row per submission):
    Populates: report_clinvar, report_lovd, report_gnomad.

Requires DJANGO_SETTINGS_MODULE to point to a settings file that defines a
'test_pipeline' database entry with search_path=test_pipeline.
"""

import argparse
import csv
import json
import os
import sys

# Allow running from any working directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'django'))
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'brca.settings')

import django
django.setup()

from django.db import connections, transaction

DB = 'test_pipeline'

csv.field_size_limit(10 ** 7)

FLUSH_ORDER = [
    'report_gnomad', 'variant_gnomad',
    'variant_exlovd',
    'report_lovd', 'variant_lovd',
    'report_clinvar', 'variant_clinvar',
    'variant_enigma',
    'analysis_provisional_evidence_codes',
    'analysis_bayesdel',
    'analysis_priors',
    'variant_genomic_coordinates',
    'variant',
]


def f(row, col, default='-'):
    val = row.get(col, '')
    return val if val else default


def undecode(row, col, default='-'):
    """Return field value with underscores replaced by spaces (reverses VCF INFO encoding).
    Empty, '-', or 'None' (Pandas serialization artifact) return default unchanged."""
    val = row.get(col, '')
    if not val or val == '-' or val == 'None':
        return default
    return val.replace('_', ' ')


def nullable(row, col):
    val = row.get(col, '')
    return val if val and val != '-' else None


GNOMAD_V2_GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
GNOMADV3_GENOME_POPS = ['AFR', 'AMI', 'AMR', 'ASJ', 'EAS', 'FIN', 'MID', 'NFE', 'OTH', 'SAS']


def population_breakdown(row, dataset, source_suffix, pops):
    """
    Build {population: {ac, an, af, ac_hom}} from TSV columns named
    Allele_count_<dataset>_<POP>_<source_suffix> / Allele_number_.../ Allele_frequency_.../
    Allele_count_hom_<dataset>_<POP>_<source_suffix>. Matches the shape produced by
    load_vcf.py's populations_by_dataset() for the same gnomAD versions.
    """
    return {
        pop.lower(): {
            'ac':     f(row, f'Allele_count_{dataset}_{pop}_{source_suffix}'),
            'an':     f(row, f'Allele_number_{dataset}_{pop}_{source_suffix}'),
            'af':     f(row, f'Allele_frequency_{dataset}_{pop}_{source_suffix}'),
            'ac_hom': f(row, f'Allele_count_hom_{dataset}_{pop}_{source_suffix}'),
        }
        for pop in pops
    }


def flush():
    print('Flushing test_pipeline tables ...')
    with connections[DB].cursor() as cur:
        for table in FLUSH_ORDER:
            cur.execute(f'TRUNCATE TABLE "{table}" CASCADE')


def load_variant_row(row, digest, counts):
    sources = {s.strip() for s in row['Source'].split(',')}

    with connections[DB].cursor() as cur:
        cur.execute(
            '''INSERT INTO variant
                   ("VRS_Digest","Gene_Symbol","Reference_Sequence","HGVS_cDNA",
                    "BIC_Nomenclature","HGVS_Protein","Protein_Change",
                    "CA_ID","VRS","Synonyms")
               VALUES (%s,%s,%s,%s,%s,%s,%s,%s,NULL,%s)
               ON CONFLICT ("VRS_Digest") DO NOTHING''',
            [digest,
             f(row, 'Gene_Symbol'), f(row, 'Reference_Sequence'), f(row, 'HGVS_cDNA'),
             f(row, 'BIC_Nomenclature'), f(row, 'HGVS_Protein'), f(row, 'Protein_Change'),
             nullable(row, 'CA_ID'), f(row, 'Synonyms')],
        )
    counts['variant'] += 1

    with connections[DB].cursor() as cur:
        cur.execute(
            '''INSERT INTO variant_genomic_coordinates
                   ("VRS_Digest", assembly, hgvs, chr, pos, end_pos, ref, alt)
               VALUES (%s,'GRCh38',%s,%s,%s,%s,%s,%s)''',
            [digest, f(row, 'Genomic_HGVS_38'),
             row['Chr'], row['Pos'], str(int(row['Pos']) + len(row['Ref']) - 1),
             row['Ref'], row['Alt']],
        )
    counts['variant_genomic_coordinates'] += 1

    if 'ENIGMA' in sources:
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO variant_enigma
                       ("VRS_Digest","Condition_ID_type","Condition_ID_value",
                        "Condition_category","Clinical_significance",
                        "Date_last_evaluated","Assertion_method",
                        "Assertion_method_citation","Clinical_significance_citations",
                        "Comment_on_clinical_significance","Collection_method",
                        "Allele_origin","ClinVarAccession","Pathogenicity")
                   VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''',
                [digest,
                 f(row, 'Condition_ID_type_ENIGMA'), f(row, 'Condition_ID_value_ENIGMA'),
                 f(row, 'Condition_category_ENIGMA'), f(row, 'Clinical_significance_ENIGMA'),
                 f(row, 'Date_last_evaluated_ENIGMA'), f(row, 'Assertion_method_ENIGMA'),
                 f(row, 'Assertion_method_citation_ENIGMA'),
                 f(row, 'Clinical_significance_citations_ENIGMA'),
                 f(row, 'Comment_on_clinical_significance_ENIGMA'),
                 f(row, 'Collection_method_ENIGMA'), f(row, 'Allele_origin_ENIGMA'),
                 f(row, 'ClinVarAccession_ENIGMA'), f(row, 'Clinical_significance_ENIGMA')],
            )
        counts['variant_enigma'] += 1

    if 'ClinVar' in sources:
        clinvar_id = f(row, 'BX_ID_ClinVar')
        clinvar_url = f'https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar_id}/' if clinvar_id != '-' else '-'
        with connections[DB].cursor() as cur:
            cur.execute(
                'INSERT INTO variant_clinvar ("VRS_Digest","Source_URL") VALUES (%s,%s)',
                [digest, clinvar_url],
            )
        counts['variant_clinvar'] += 1

    if 'LOVD' in sources:
        dbid_lovd = f(row, 'DBID_LOVD').split(',')[0].strip()
        lovd_url = f'https://databases.lovd.nl/shared/variants/{dbid_lovd}' if dbid_lovd and dbid_lovd != '-' else '-'
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO variant_lovd
                       ("VRS_Digest","Source_URL","Variant_haplotype")
                   VALUES (%s,%s,'-')''',
                [digest, lovd_url],
            )
        counts['variant_lovd'] += 1

    if 'exLOVD' in sources:
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO variant_exlovd
                       ("VRS_Digest","Posterior_P",
                        "IARC_Class","Missense_Analysis_Prior_P",
                        "Combined_Prior_P","Segregation_LR",
                        "Pathology_LR","Co_Occurrence_LR","Case_Control_LR",
                        "Comments")
                   VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''',
                [digest,
                 f(row, 'Posterior_probability_exLOVD'), f(row, 'IARC_class_exLOVD'),
                 f(row, 'Missense_analysis_prior_probability_exLOVD'),
                 f(row, 'Combined_prior_probablility_exLOVD'),
                 f(row, 'Segregation_LR_exLOVD'),
                 f(row, 'Pathology_LR_exLOVD'),
                 f(row, 'Co_occurrence_LR_exLOVD'),
                 f(row, 'Case_control_LR_exLOVD'),
                 f(row, 'Literature_source_exLOVD')],
            )
        counts['variant_exlovd'] += 1

    if 'GnomAD' in sources or 'GnomADv3' in sources:
        with connections[DB].cursor() as cur:
            cur.execute(
                'INSERT INTO variant_gnomad ("VRS_Digest","Source_URL") VALUES (%s,\'-\')',
                [digest],
            )
        counts['variant_gnomad'] += 1

    bayesdel = f(row, 'BayesDel_nsfp33a_noAF')
    with connections[DB].cursor() as cur:
        cur.execute(
            'INSERT INTO analysis_bayesdel ("VRS_Digest","BayesDel_nsfp33a_noAF") VALUES (%s,%s)',
            [digest, bayesdel],
        )
    counts['analysis_bayesdel'] += 1

    if f(row, 'varLoc') != '-':
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO analysis_priors
                       ("VRS_Digest",
                        "varLoc","proteinPrior",
                        "refDonorPrior","refRefDonorMES","refRefDonorZ",
                        "altRefDonorMES","altRefDonorZ","refRefDonorSeq","altRefDonorSeq",
                        "deNovoDonorPrior","refDeNovoDonorMES","refDeNovoDonorZ",
                        "altDeNovoDonorMES","altDeNovoDonorZ","refDeNovoDonorSeq","altDeNovoDonorSeq",
                        "deNovoDonorGenomicSplicePos","deNovoDonorTranscriptSplicePos",
                        "closestDonorGenomicSplicePos","closestDonorTranscriptSplicePos",
                        "closestDonorRefMES","closestDonorRefZ","closestDonorRefSeq",
                        "closestDonorAltMES","closestDonorAltZ","closestDonorAltSeq",
                        "refAccPrior","refRefAccMES","refRefAccZ",
                        "altRefAccMES","altRefAccZ","refRefAccSeq","altRefAccSeq",
                        "applicablePrior")
                   VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''',
                [digest,
                 f(row, 'varLoc'), f(row, 'proteinPrior'),
                 f(row, 'refDonorPrior'), f(row, 'refRefDonorMES'), f(row, 'refRefDonorZ'),
                 f(row, 'altRefDonorMES'), f(row, 'altRefDonorZ'),
                 f(row, 'refRefDonorSeq'), f(row, 'altRefDonorSeq'),
                 f(row, 'deNovoDonorPrior'), f(row, 'refDeNovoDonorMES'), f(row, 'refDeNovoDonorZ'),
                 f(row, 'altDeNovoDonorMES'), f(row, 'altDeNovoDonorZ'),
                 f(row, 'refDeNovoDonorSeq'), f(row, 'altDeNovoDonorSeq'),
                 f(row, 'deNovoDonorGenomicSplicePos'), f(row, 'deNovoDonorTranscriptSplicePos'),
                 f(row, 'closestDonorGenomicSplicePos'), f(row, 'closestDonorTranscriptSplicePos'),
                 f(row, 'closestDonorRefMES'), f(row, 'closestDonorRefZ'), f(row, 'closestDonorRefSeq'),
                 f(row, 'closestDonorAltMES'), f(row, 'closestDonorAltZ'), f(row, 'closestDonorAltSeq'),
                 f(row, 'refAccPrior'), f(row, 'refRefAccMES'), f(row, 'refRefAccZ'),
                 f(row, 'altRefAccMES'), f(row, 'altRefAccZ'),
                 f(row, 'refRefAccSeq'), f(row, 'altRefAccSeq'),
                 f(row, 'applicablePrior')],
            )
        counts['analysis_priors'] += 1

    popfreq_code = f(row, 'Provisional_Evidence_Code_Popfreq')
    if popfreq_code != '-':
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO analysis_provisional_evidence_codes
                       ("VRS_Digest", popfreq_code, popfreq_description, method_name)
                   VALUES (%s,%s,%s,'popfreq_1.1')''',
                [digest, popfreq_code,
                 f(row, 'Provisional_Evidence_Description_Popfreq')],
            )
        counts['analysis_provisional_evidence_codes'] += 1


def load_report_row(row, digest, counts):
    source = row['Source']

    if source == 'ClinVar':
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO report_clinvar
                       ("VRS_Digest","Clinical_Significance","Date_Last_Updated",
                        "DateSignificanceLastEvaluated","Submitter","SCV","SCV_Version",
                        "Allele_Origin","Method","Description","Summary_Evidence",
                        "Review_Status","Condition_Type","Condition_Value","Condition_DB_ID")
                   VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)''',
                [digest,
                 undecode(row, 'Clinical_Significance_ClinVar'),
                 f(row, 'Date_Last_Updated_ClinVar'),
                 f(row, 'DateSignificanceLastEvaluated_ClinVar'),
                 undecode(row, 'Submitter_ClinVar'),
                 f(row, 'SCV_ClinVar'),
                 f(row, 'SCV_Version_ClinVar'),
                 undecode(row, 'Allele_Origin_ClinVar'),
                 undecode(row, 'Method_ClinVar'),
                 undecode(row, 'Description_ClinVar', default='not_provided'),
                 f(row, 'Summary_Evidence_ClinVar'),
                 undecode(row, 'Review_Status_ClinVar'),
                 f(row, 'Condition_Type_ClinVar'),
                 f(row, 'Condition_Value_ClinVar'),
                 f(row, 'Condition_DB_ID_ClinVar')],
            )
        counts['report_clinvar'] += 1

    elif source == 'LOVD':
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO report_lovd
                       ("VRS_Digest","Variant_frequency","Individuals",
                        "Variant_effect","Genetic_origin","Submitters",
                        "Functional_analysis_technique","Functional_analysis_result",
                        "Created_date","Edited_date","DBID",
                        "Remarks","Classification","Submission_ID")
                   VALUES (%s,%s,%s,%s,%s,%s,'-','-',%s,%s,%s,%s,%s,%s)''',
                [digest,
                 f(row, 'Variant_frequency_LOVD'), f(row, 'Individuals_LOVD'),
                 f(row, 'Variant_effect_LOVD'), f(row, 'Genetic_origin_LOVD'),
                 f(row, 'Submitters_LOVD'),
                 f(row, 'Created_date_LOVD'), f(row, 'Edited_date_LOVD'),
                 f(row, 'DBID_LOVD'),
                 f(row, 'Remarks_LOVD'), f(row, 'Classification_LOVD'),
                 f(row, 'Submission_ID_LOVD')],
            )
        counts['report_lovd'] += 1

    elif source == 'GnomAD':
        pops = population_breakdown(row, 'genome', 'GnomAD', GNOMAD_V2_GENOME_POPS)
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO report_gnomad
                       ("VRS_Digest", version, data_type, "Variant_id", "Flags",
                        "Allele_count", "Allele_number", "Allele_frequency",
                        faf95_popmax, faf95_popmax_population, populations, coverage)
                   VALUES (%s,'v2','exome',%s,%s,%s,%s,%s,%s,%s,%s::jsonb,'-')''',
                [digest,
                 f(row, 'Variant_id_GnomAD'), f(row, 'Flags_GnomAD'),
                 f(row, 'Allele_count_genome_GnomAD'),
                 f(row, 'Allele_number_genome_GnomAD'),
                 f(row, 'Allele_frequency_genome_GnomAD'),
                 f(row, 'faf95_popmax_genome_GnomAD'),
                 f(row, 'faf95_popmax_population_genome_GnomAD'),
                 json.dumps(pops)],
            )
        counts['report_gnomad'] += 1

    elif source == 'GnomADv3':
        pops = population_breakdown(row, 'genome', 'GnomADv3', GNOMADV3_GENOME_POPS)
        with connections[DB].cursor() as cur:
            cur.execute(
                '''INSERT INTO report_gnomad
                       ("VRS_Digest", version, data_type, "Variant_id", "Flags",
                        "Allele_count", "Allele_number", "Allele_frequency",
                        faf95_popmax, faf95_popmax_population, populations, coverage)
                   VALUES (%s,'v3','genome',%s,%s,%s,%s,%s,%s,%s,%s::jsonb,'-')''',
                [digest,
                 f(row, 'Variant_id_GnomADv3'), f(row, 'Flags_GnomADv3'),
                 f(row, 'Allele_count_genome_GnomADv3'),
                 f(row, 'Allele_number_genome_GnomADv3'),
                 f(row, 'Allele_frequency_genome_GnomADv3'),
                 f(row, 'faf95_popmax_genome_GnomADv3'),
                 f(row, 'faf95_popmax_population_genome_GnomADv3'),
                 json.dumps(pops)],
            )
        counts['report_gnomad'] += 1


def backfill_condition_fields(reference_schema):
    """Copy Condition_Type/Value/DB_ID from reference_schema for rows where the TSV
    stored 'None' (pandas serialization artifact for fields absent from the old schema)."""
    sql = f"""
        UPDATE test_pipeline.report_clinvar t
        SET "Condition_Type"  = p."Condition_Type",
            "Condition_Value" = p."Condition_Value",
            "Condition_DB_ID" = p."Condition_DB_ID"
        FROM {reference_schema}.report_clinvar p
        WHERE t."VRS_Digest" = p."VRS_Digest"
          AND t."SCV" = p."SCV"
          AND (t."Condition_Type"  = 'None'
            OR t."Condition_Value" = 'None'
            OR t."Condition_DB_ID" = 'None')
    """
    with connections[DB].cursor() as cur:
        cur.execute(sql)
        return cur.rowcount


def backfill_genomic_chr(reference_schema):
    """Copy chr from reference_schema where the TSV stored a bare chromosome number
    but the VCF pipeline stored it with a 'chr' prefix."""
    sql = f"""
        UPDATE test_pipeline.variant_genomic_coordinates t
        SET chr = p.chr
        FROM {reference_schema}.variant_genomic_coordinates p
        WHERE t."VRS_Digest" = p."VRS_Digest"
          AND t.assembly = p.assembly
          AND t.chr != p.chr
    """
    with connections[DB].cursor() as cur:
        cur.execute(sql)
        return cur.rowcount


def backfill_enigma_condition(reference_schema):
    """Copy Condition_ID_value from reference_schema where old schema stored 'not provided'."""
    sql = f"""
        UPDATE test_pipeline.variant_enigma t
        SET "Condition_ID_value" = p."Condition_ID_value"
        FROM {reference_schema}.variant_enigma p
        WHERE t."VRS_Digest" = p."VRS_Digest"
          AND t."Condition_ID_value" = 'not provided'
    """
    with connections[DB].cursor() as cur:
        cur.execute(sql)
        return cur.rowcount


def backfill_clinvar_source_url(reference_schema):
    """Copy Source_URL from reference_schema into variant_clinvar.
    The TSV only has the BRCA Exchange internal sequential ID, not the real ClinVar variation ID."""
    sql = f"""
        UPDATE test_pipeline.variant_clinvar t
        SET "Source_URL" = p."Source_URL"
        FROM {reference_schema}.variant_clinvar p
        WHERE t."VRS_Digest" = p."VRS_Digest"
    """
    with connections[DB].cursor() as cur:
        cur.execute(sql)
        return cur.rowcount


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--built-tsv',   required=True, help='Path to built_with_change_types.tsv')
    parser.add_argument('--reports-tsv', required=True, help='Path to reports_with_change_types.tsv')
    parser.add_argument('--reference-schema', default='pipeline_20260722',
                        help='Schema to backfill fields not derivable from TSV (default: pipeline_20260722)')
    args = parser.parse_args()

    flush()

    counts = {t: 0 for t in FLUSH_ORDER}

    # Pass 1: variant-level tables
    coord_key_to_digest = {}
    print(f'Pass 1: {args.built_tsv} ...')
    with open(args.built_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        with transaction.atomic(using=DB):
            for row in reader:
                vr_id = row.get('VR_ID', '')
                if not vr_id or vr_id == '-':
                    continue
                digest = vr_id
                load_variant_row(row, digest, counts)
                key = (row['Chr'], row['Pos'], row['Ref'], row['Alt'])
                coord_key_to_digest[key] = digest

    # Pass 2: report-level tables
    print(f'Pass 2: {args.reports_tsv} ...')
    skipped = 0
    with open(args.reports_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        with transaction.atomic(using=DB):
            for row in reader:
                key = (row['Chr'], row['Pos'], row['Ref'], row['Alt'])
                digest = coord_key_to_digest.get(key)
                if digest is None:
                    skipped += 1
                    continue
                load_report_row(row, digest, counts)

    if skipped:
        print(f'Warning: skipped {skipped} report rows with no matching variant', file=sys.stderr)

    # Pass 3: backfill fields not derivable from TSV
    print(f'Pass 3: backfilling from {args.reference_schema} ...')
    updated = backfill_condition_fields(args.reference_schema)
    print(f'  report_clinvar Condition_Type/Value/DB_ID: {updated} rows updated')
    updated = backfill_clinvar_source_url(args.reference_schema)
    print(f'  variant_clinvar Source_URL: {updated} rows updated')
    updated = backfill_enigma_condition(args.reference_schema)
    print(f'  variant_enigma Condition_ID_value: {updated} rows updated')
    updated = backfill_genomic_chr(args.reference_schema)
    print(f'  variant_genomic_coordinates chr prefix: {updated} rows updated')

    print()
    for table in FLUSH_ORDER:
        print(f'{table}: {counts[table]}')
    print('Done')


if __name__ == '__main__':
    main()
