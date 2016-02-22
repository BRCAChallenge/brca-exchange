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

        SELECT regexp_split_to_table("Gene_symbol", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Genomic_Coordinate", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_genomic", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Transcript_id", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_cDNA", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_protein", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("BIC_Nomenclature", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Abbrev_AA_change", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("URL", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Condition_ID_type", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Condition_ID_value", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Condition_category", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Clinical_significance", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Date_last_evaluated", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Assertion_method", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Assertion_method_citation", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Clinical_significance_citations", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Comment_on_clinical_significance", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Collection_method", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Allele_origin", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("ClinVarAccession", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("SAS_Allele_frequency_1000_Genomes", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("EAS_Allele_frequency_1000_Genomes", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Allele_frequency_1000_Genomes", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("AMR_Allele_frequency_1000Genomes", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("EUR_Allele_frequency_1000Genomes", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("AFR_Allele_frequency_1000Genomes", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Allele_origin_ClinVar", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Variant_clinical_significance_ClinVar", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_genomic_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Origin_of_variant_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_protein_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Variant_frequency_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_cDNA_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Variant_affecting_protein_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Variant_haplotype_LOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("VEP_Gene_ExAC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Allele_frequency_ExAC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("VEP_HGVSc_ExAC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("VEP_Consequence_ExAC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("VEP_HGVSp_ExAC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Exon_number_exLOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("IARC_class_exLOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("BIC_exLOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_cDNA_exLOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Literature_source_exLOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_protein_exLOVD", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Number_of_family_member_carrying_mutation_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_genomic_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Germline_or_Somatic_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Mutation_type_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("BIC_Designation_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Literature_citation_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Exon_number_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Clinical_importance_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Clinical_classification_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_protein_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("HGVS_cDNA_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Ethnicity_BIC", '[\s+|]') as word from variant UNION
        SELECT regexp_split_to_table("Patient_nationality_BIC", '[\s+|]') as word from variant

        ) AS combined_words;

        CREATE INDEX words_idx ON words USING gin(word gin_trgm_ops);
    """)

    ]
