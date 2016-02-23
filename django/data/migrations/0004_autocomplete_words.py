from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0003_data'),
    ]

    operations = [
        migrations.RunSQL("""
        CREATE EXTENSION pg_trgm;

        CREATE TABLE words AS SELECT DISTINCT word FROM (

        SELECT regexp_split_to_table(lower("Gene_symbol"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Genomic_Coordinate"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_genomic"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Transcript_id"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_cDNA"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_protein"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("BIC_Nomenclature"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Abbrev_AA_change"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Condition_ID_type"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Condition_ID_value"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Condition_category"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Clinical_significance"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Assertion_method"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Clinical_significance_citations"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Allele_origin"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("ClinVarAccession"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_genomic_LOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Origin_of_variant_LOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_protein_LOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_cDNA_LOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Variant_haplotype_LOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("VEP_Gene_ExAC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("VEP_HGVSc_ExAC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("VEP_Consequence_ExAC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("VEP_HGVSp_ExAC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("IARC_class_exLOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("BIC_exLOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_cDNA_exLOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_protein_exLOVD"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_genomic_BIC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Mutation_type_BIC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("BIC_Designation_BIC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_protein_BIC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("HGVS_cDNA_BIC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Ethnicity_BIC"), '[\s+|''"]') as word from variant UNION
        SELECT regexp_split_to_table(lower("Patient_nationality_BIC"), '[\s+|''"]') as word from variant

        ) AS combined_words;

        CREATE INDEX words_idx ON words USING gin(word gin_trgm_ops);
    """)

    ]
