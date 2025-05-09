# -*- coding: utf-8 -*-
# Generated by Django 1.9.13 on 2025-04-22 15:43
from __future__ import unicode_literals

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0058_add_provisional_evidence_fields_to_variant'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='variant',
            name='Discordant',
        ),
    
        migrations.RunSQL("""
            ALTER TABLE variant DROP COLUMN fts_document CASCADE;
            ALTER TABLE variant DROP COLUMN fts_standard CASCADE;
            DROP FUNCTION variant_fts_standard(v variant) CASCADE;
            DROP FUNCTION variant_fts_document(v variant) CASCADE;
            DROP FUNCTION variant_fts_trigger() CASCADE;

            ALTER TABLE variant ADD COLUMN fts_document TSVECTOR;
            ALTER TABLE variant ADD COLUMN fts_standard TSVECTOR;

            CREATE FUNCTION variant_fts_standard(v variant) RETURNS tsvector AS $$
                DECLARE
                    fts_standard TEXT;
                BEGIN
                    SELECT concat_ws(' ',
                        v."Source",
                        v."URL_ENIGMA",
                        v."Condition_ID_type_ENIGMA",
                        v."Condition_ID_value_ENIGMA",
                        v."Condition_category_ENIGMA",
                        v."Clinical_significance_ENIGMA",
                        v."Date_last_evaluated_ENIGMA",
                        v."Assertion_method_ENIGMA",
                        v."Assertion_method_citation_ENIGMA",
                        v."Clinical_significance_citations_ENIGMA",
                        v."Comment_on_clinical_significance_ENIGMA",
                        v."Collection_method_ENIGMA",
                        v."Allele_origin_ENIGMA",
                        v."ClinVarAccession_ENIGMA",
                        v."Clinical_Significance_ClinVar",
                        v."Date_Last_Updated_ClinVar",
                        v."Submitter_ClinVar",
                        v."SCV_ClinVar",
                        v."Allele_Origin_ClinVar",
                        v."Method_ClinVar",
                        v."Functional_analysis_result_LOVD",
                        v."Functional_analysis_technique_LOVD",
                        v."Variant_frequency_LOVD",
                        v."Variant_haplotype_LOVD",
                        v."Minor_allele_frequency_percent_ESP",
                        v."EUR_Allele_frequency_1000_Genomes",
                        v."AFR_Allele_frequency_1000_Genomes",
                        v."AMR_Allele_frequency_1000_Genomes",
                        v."EAS_Allele_frequency_1000_Genomes",
                        v."Allele_frequency_1000_Genomes",
                        v."SAS_Allele_frequency_1000_Genomes",
                        v."Allele_frequency_ExAC",
                        v."Patient_nationality_BIC",
                        v."Clinical_importance_BIC",
                        v."Clinical_classification_BIC",
                        v."Literature_citation_BIC",
                        v."Number_of_family_member_carrying_mutation_BIC",
                        v."Germline_or_Somatic_BIC",
                        v."Ethnicity_BIC",
                        v."Mutation_type_BIC",
                        v."IARC_class_exLOVD",
                        v."Sum_family_LR_exLOVD",
                        v."Combined_prior_probablility_exLOVD",
                        v."Literature_source_exLOVD",
                        v."Co_occurrence_LR_exLOVD",
                        v."Posterior_probability_exLOVD",
                        v."Missense_analysis_prior_probability_exLOVD",
                        v."Segregation_LR_exLOVD",
                        v."Gene_Symbol",
                        v."Reference_Sequence",
                        v."HGVS_cDNA",
                        v."BIC_Nomenclature",
                        v."HGVS_Protein",
                        v."HGVS_RNA",
                        v."Protein_Change",
                        v."Allele_Frequency",
                        v."Max_Allele_Frequency",
                        v."Genomic_Coordinate_hg38",
                        v."Source_URL",
                        v."Pathogenicity_expert",
                        v."Pathogenicity_all")
            INTO fts_standard;
                    RETURN to_tsvector('pg_catalog.simple', fts_standard);
            END;
            $$ LANGUAGE plpgsql;

            CREATE FUNCTION variant_fts_document(v variant) RETURNS tsvector AS $$
                BEGIN
                    RETURN variant_fts_standard(v) || variant_fts_synonyms(v);
                END;
                $$ LANGUAGE plpgsql;

            CREATE FUNCTION variant_fts_trigger() RETURNS TRIGGER AS $$
            BEGIN
            NEW.fts_standard=variant_fts_standard(NEW);
            NEW.fts_document=variant_fts_document(NEW);
            RETURN NEW;
            END;
            $$ LANGUAGE plpgsql;

            CREATE TRIGGER variant_fts_update_trigger BEFORE UPDATE ON variant FOR EACH ROW EXECUTE PROCEDURE variant_fts_trigger();
            CREATE TRIGGER variant_fts_insert_trigger BEFORE INSERT ON variant FOR EACH ROW EXECUTE PROCEDURE variant_fts_trigger();

            CREATE INDEX variant_fts_document_index ON variant USING gin(fts_document);
            CREATE INDEX variant_fts_standard_index ON variant USING gin(fts_standard);

            """),
        migrations.RunSQL(
            """
            DROP MATERIALIZED VIEW IF EXISTS currentvariant;
            CREATE MATERIALIZED VIEW currentvariant AS (
                SELECT * FROM "variant" WHERE (
                    "id" IN ( SELECT DISTINCT ON ("Genomic_Coordinate_hg38") "id" FROM "variant" ORDER BY "Genomic_Coordinate_hg38" ASC, "Data_Release_id" DESC )
                )
            );
            """
        )
    ]
