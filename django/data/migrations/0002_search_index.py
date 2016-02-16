from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):
    dependencies = [
        ('data', '0001_initial'),
    ]

    operations = [
        migrations.RunSQL("""
    ALTER TABLE variant ADD COLUMN fts_document TSVECTOR;

    CREATE FUNCTION variant_fts_document(v variant) RETURNS tsvector AS $$
        DECLARE
            fts_document TEXT;
        BEGIN
            SELECT concat_ws(' ',
                v."Gene_symbol",
                v."Genomic_Coordinate",
                v."HGVS_genomic",
                v."Transcript_id",
                v."HGVS_cDNA",
                v."HGVS_protein",
                v."BIC_Nomenclature",
                v."Abbrev_AA_change",
                v."URL",
                v."Condition_ID_type",
                v."Condition_ID_value",
                v."Condition_category",
                v."Clinical_significance",
                v."Date_last_evaluated",
                v."Assertion_method",
                v."Assertion_method_citation",
                v."Clinical_significance_citations",
                v."Comment_on_clinical_significance",
                v."Collection_method",
                v."Allele_origin",
                v."ClinVarAccession",
                v."SAS_Allele_frequency_1000_Genomes",
                v."EAS_Allele_frequency_1000_Genomes",
                v."Allele_frequency_1000_Genomes",
                v."AMR_Allele_frequency_1000Genomes",
                v."EUR_Allele_frequency_1000Genomes",
                v."AFR_Allele_frequency_1000Genomes",
                v."Allele_origin_ClinVar",
                v."Variant_clinical_significance_ClinVar",
                v."HGVS_genomic_LOVD",
                v."Origin_of_variant_LOVD",
                v."HGVS_protein_LOVD",
                v."Variant_frequency_LOVD",
                v."HGVS_cDNA_LOVD",
                v."Variant_affecting_protein_LOVD",
                v."Variant_haplotype_LOVD",
                v."VEP_Gene_ExAC",
                v."Allele_frequency_ExAC",
                v."VEP_HGVSc_ExAC",
                v."VEP_Consequence_ExAC",
                v."VEP_HGVSp_ExAC",
                v."Exon_number_exLOVD",
                v."IARC_class_exLOVD",
                v."BIC_exLOVD",
                v."HGVS_cDNA_exLOVD",
                v."Literature_source_exLOVD",
                v."HGVS_protein_exLOVD",
                v."Number_of_family_member_carrying_mutation_BIC",
                v."HGVS_genomic_BIC",
                v."Germline_or_Somatic_BIC",
                v."Mutation_type_BIC",
                v."BIC_Designation_BIC",
                v."Literature_citation_BIC",
                v."Exon_number_BIC",
                v."Clinical_importance_BIC",
                v."Clinical_classification_BIC",
                v."HGVS_protein_BIC",
                v."HGVS_cDNA_BIC",
                v."Ethnicity_BIC",
                v."Patient_nationality_BIC")
    INTO fts_document;
            RETURN to_tsvector('pg_catalog.simple', fts_document);
    END;
    $$ LANGUAGE plpgsql;

    CREATE FUNCTION variant_fts_document_trigger() RETURNS TRIGGER AS $$
    BEGIN
    NEW.fts_document=variant_fts_document(NEW);
    RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;

    CREATE TRIGGER variant_fts_update_trigger BEFORE UPDATE ON variant FOR EACH ROW EXECUTE PROCEDURE variant_fts_document_trigger();
    CREATE TRIGGER variant_fts_insert_trigger BEFORE INSERT ON variant FOR EACH ROW EXECUTE PROCEDURE variant_fts_document_trigger();

    CREATE INDEX variant_fts_index ON variant USING gin(fts_document);

    """)

    ]
