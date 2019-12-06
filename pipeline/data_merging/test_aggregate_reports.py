import unittest
import os
from os import path
import aggregate_reports
import utilities

VCF_TESTDATA_FILENAME = path.join(path.dirname(__file__), 'test_files/1000_Genomes.vcf')
TSV_TESTDATA_FILENAME = path.join(path.dirname(__file__), 'test_files/enigma_from_clinvar.tsv')
INPUT_DIRECTORY = path.join(path.dirname(__file__), 'test_files/')


class TestStringMethods(unittest.TestCase):

    def setUp(self):
        self.columns = ['Source', 'Gene_symbol_ENIGMA', 'Genomic_Coordinate',
                        'Chr', 'Pos', 'Ref', 'Alt', 'Reference_sequence_ENIGMA',
                        'HGVS_cDNA_ENIGMA', 'BIC_Nomenclature_ENIGMA',
                        'Abbrev_AA_change_ENIGMA', 'URL_ENIGMA',
                        'Condition_ID_type_ENIGMA', 'Condition_ID_value_ENIGMA',
                        'Condition_category_ENIGMA',
                        'Clinical_significance_ENIGMA',
                        'Date_last_evaluated_ENIGMA', 'Assertion_method_ENIGMA',
                        'Assertion_method_citation_ENIGMA',
                        'Clinical_significance_citations_ENIGMA',
                        'Comment_on_clinical_significance_ENIGMA',
                        'Collection_method_ENIGMA', 'Allele_origin_ENIGMA',
                        'ClinVarAccession_ENIGMA', 'HGVS_protein_ENIGMA',
                        'BX_ID_ENIGMA', 'Clinical_Significance_ClinVar',
                        'Date_Last_Updated_ClinVar', 'BX_ID_ClinVar',
                        'HGVS_ClinVar', 'Submitter_ClinVar', 'Protein_ClinVar',
                        'SCV_ClinVar', 'Allele_Origin_ClinVar',
                        'Method_ClinVar', 'Description_ClinVar',
                        'Summary_Evidence_ClinVar', 'Review_Status_ClinVar',
                        'Condition_Type_ClinVar', 'Condition_Value_ClinVar',
                        'Condition_DB_ID_ClinVar', 'Individuals_LOVD',
                        'BX_ID_LOVD', 'Variant_effect_LOVD', 'Variant_frequency_LOVD',
                        'HGVS_cDNA_LOVD', 'HGVS_protein_LOVD', 'Genetic_origin_LOVD',
                        'RNA_LOVD', 'Submitters_LOVD', 'DBID_LOVD',
                        'Created_date_LOVD',
                        'Edited_date_LOVD', 'Submission_ID_LOVD',
                        'Remarks_LOVD', 'Classification_LOVD',
                        'BX_ID_ESP', 'Minor_allele_frequency_percent_ESP',
                        'EA_Allele_Frequency_ESP', 'AA_Allele_Frequency_ESP',
                        'Allele_Frequency_ESP', 'polyPhen2_result_ESP',
                        'EUR_Allele_frequency_1000_Genomes',
                        'AFR_Allele_frequency_1000_Genomes',
                        'AMR_Allele_frequency_1000_Genomes',
                        'EAS_Allele_frequency_1000_Genomes',
                        'BX_ID_1000_Genomes', 'Allele_frequency_1000_Genomes',
                        'SAS_Allele_frequency_1000_Genomes',
                        'Allele_frequency_ExAC', 'BX_ID_ExAC', 'BX_ID_BIC',
                        'Patient_nationality_BIC', 'Clinical_importance_BIC',
                        'Clinical_classification_BIC', 'BIC_Designation_BIC',
                        'Literature_citation_BIC',
                        'Number_of_family_member_carrying_mutation_BIC',
                        'Germline_or_Somatic_BIC', 'Ethnicity_BIC',
                        'Mutation_type_BIC', 'IARC_class_exLOVD',
                        'BIC_Nomenclature_exLOVD', 'Sum_family_LR_exLOVD',
                        'Combined_prior_probablility_exLOVD', 'BX_ID_exLOVD',
                        'HGVS_cDNA_exLOVD', 'Literature_source_exLOVD',
                        'Co_occurrence_LR_exLOVD',
                        'Posterior_probability_exLOVD',
                        'Missense_analysis_prior_probability_exLOVD',
                        'Segregation_LR_exLOVD', 'HGVS_protein_exLOVD',
                        "Allele_count_AFR_ExAC", "Allele_number_AFR_ExAC",
                        "Homozygous_count_AFR_ExAC", "Allele_count_AMR_ExAC",
                        "Allele_number_AMR_ExAC", "Homozygous_count_AMR_ExAC",
                        "Allele_count_EAS_ExAC", "Allele_number_EAS_ExAC",
                        "Homozygous_count_EAS_ExAC", "Allele_count_FIN_ExAC",
                        "Allele_number_FIN_ExAC", "Homozygous_count_FIN_ExAC",
                        "Allele_count_NFE_ExAC", "Allele_number_NFE_ExAC",
                        "Homozygous_count_NFE_ExAC", "Allele_count_OTH_ExAC",
                        "Allele_number_OTH_ExAC", "Homozygous_count_OTH_ExAC",
                        "Allele_count_SAS_ExAC", "Allele_number_SAS_ExAC",
                        "Homozygous_count_SAS_ExAC",
                        "Allele_frequency_AFR_ExAC",
                        "Allele_frequency_AMR_ExAC",
                        "Allele_frequency_EAS_ExAC",
                        "Allele_frequency_FIN_ExAC",
                        "Allele_frequency_NFE_ExAC",
                        "Allele_frequency_OTH_ExAC",
                        "Allele_frequency_SAS_ExAC",
                        "DateSignificanceLastEvaluated_ClinVar",
                        "SCV_Version_ClinVar",
                        "Synonyms_ClinVar",
                        "HGVS_Nucleotide_Findlay_BRCA1_Ring_Function_Scores",
                        "Log_RNA_Depletion_Findlay_BRCA1_Ring_Function_Scores",
                        "Functional_Enrichment_Score_Findlay_BRCA1_Ring_Function_Scores",
                        "BX_ID_Findlay_BRCA1_Ring_Function_Scores",
                        "HGVS_cDNA_GnomAD",
                        "HGVS_GnomAD",
                        "HGVS_protein_GnomAD",
                        "Flags_GnomAD",
                        "Consequence_GnomAD",
                        "Variant_id_GnomAD",
                        "Allele_count_genome_AFR_GnomAD",
                        "Allele_count_hemi_genome_AFR_GnomAD",
                        "Allele_count_hom_genome_AFR_GnomAD",
                        "Allele_number_genome_AFR_GnomAD",
                        "Allele_frequency_genome_AFR_GnomAD",
                        "Allele_count_genome_AMR_GnomAD",
                        "Allele_count_hemi_genome_AMR_GnomAD",
                        "Allele_count_hom_genome_AMR_GnomAD",
                        "Allele_number_genome_AMR_GnomAD",
                        "Allele_frequency_genome_AMR_GnomAD",
                        "Allele_count_genome_ASJ_GnomAD",
                        "Allele_count_hemi_genome_ASJ_GnomAD",
                        "Allele_count_hom_genome_ASJ_GnomAD",
                        "Allele_number_genome_ASJ_GnomAD",
                        "Allele_frequency_genome_ASJ_GnomAD",
                        "Allele_count_genome_EAS_GnomAD",
                        "Allele_count_hemi_genome_EAS_GnomAD",
                        "Allele_count_hom_genome_EAS_GnomAD",
                        "Allele_number_genome_EAS_GnomAD",
                        "Allele_frequency_genome_EAS_GnomAD",
                        "Allele_count_genome_FIN_GnomAD",
                        "Allele_count_hemi_genome_FIN_GnomAD",
                        "Allele_count_hom_genome_FIN_GnomAD",
                        "Allele_number_genome_FIN_GnomAD",
                        "Allele_frequency_genome_FIN_GnomAD",
                        "Allele_count_genome_NFE_GnomAD",
                        "Allele_count_hemi_genome_NFE_GnomAD",
                        "Allele_count_hom_genome_NFE_GnomAD",
                        "Allele_number_genome_NFE_GnomAD",
                        "Allele_frequency_genome_NFE_GnomAD",
                        "Allele_count_genome_OTH_GnomAD",
                        "Allele_count_hemi_genome_OTH_GnomAD",
                        "Allele_count_hom_genome_OTH_GnomAD",
                        "Allele_number_genome_OTH_GnomAD",
                        "Allele_frequency_genome_OTH_GnomAD",
                        "Allele_count_genome_SAS_GnomAD",
                        "Allele_count_hemi_genome_SAS_GnomAD",
                        "Allele_count_hom_genome_SAS_GnomAD",
                        "Allele_number_genome_SAS_GnomAD",
                        "Allele_frequency_genome_SAS_GnomAD",
                        "Allele_count_genome_GnomAD",
                        "Allele_number_genome_GnomAD",
                        "Allele_frequency_genome_GnomAD",
                        "Allele_count_exome_AFR_GnomAD",
                        "Allele_count_exome_AFR_GnomAD",
                        "Allele_count_hom_exome_AFR_GnomAD",
                        "Allele_number_exome_AFR_GnomAD",
                        "Allele_frequency_exome_AFR_GnomAD",
                        "Allele_count_exome_AMR_GnomAD",
                        "Allele_count_hemi_exome_AMR_GnomAD",
                        "Allele_count_hom_exome_AMR_GnomAD",
                        "Allele_number_exome_AMR_GnomAD",
                        "Allele_frequency_exome_AMR_GnomAD",
                        "Allele_count_exome_ASJ_GnomAD",
                        "Allele_count_hemi_exome_ASJ_GnomAD",
                        "Allele_count_hom_exome_ASJ_GnomAD",
                        "Allele_number_exome_ASJ_GnomAD",
                        "Allele_frequency_exome_ASJ_GnomAD",
                        "Allele_count_exome_EAS_GnomAD",
                        "Allele_count_hemi_exome_EAS_GnomAD",
                        "Allele_count_hom_exome_EAS_GnomAD",
                        "Allele_number_exome_EAS_GnomAD",
                        "Allele_frequency_exome_EAS_GnomAD",
                        "Allele_count_exome_FIN_GnomAD",
                        "Allele_count_hemi_exome_FIN_GnomAD",
                        "Allele_count_hom_exome_FIN_GnomAD",
                        "Allele_number_exome_FIN_GnomAD",
                        "Allele_frequency_exome_FIN_GnomAD",
                        "Allele_count_exome_NFE_GnomAD",
                        "Allele_count_hemi_exome_NFE_GnomAD",
                        "Allele_count_hom_exome_NFE_GnomAD",
                        "Allele_number_exome_NFE_GnomAD",
                        "Allele_frequency_exome_NFE_GnomAD",
                        "Allele_count_exome_OTH_GnomAD",
                        "Allele_count_hemi_exome_OTH_GnomAD",
                        "Allele_count_hom_exome_OTH_GnomAD",
                        "Allele_number_exome_OTH_GnomAD",
                        "Allele_frequency_exome_OTH_GnomAD",
                        "Allele_count_exome_SAS_GnomAD",
                        "Allele_count_hemi_exome_SAS_GnomAD",
                        "Allele_count_hom_exome_SAS_GnomAD",
                        "Allele_number_exome_SAS_GnomAD",
                        "Allele_frequency_exome_SAS_GnomAD",
                        "Allele_number_exome_GnomAD",
                        "Allele_count_exome_GnomAD",
                        "Allele_frequency_exome_GnomAD",
                        "BX_ID_GnomAD"]

        self.sources = aggregate_reports.FIELD_DICT.keys() + ["ENIGMA"]
        self.vcf_test_file = VCF_TESTDATA_FILENAME
        self.tsv_test_file = TSV_TESTDATA_FILENAME

        pwd = os.path.dirname(os.path.realpath(__file__))

        gene_config_df = utilities.load_config(os.path.join(pwd, 'test_files', 'gene_config_test.txt'))
        self.genome_regions_symbol_dict = utilities.get_genome_regions_symbol_dict(gene_config_df, 'start_hg38_legacy_variants', 'end_hg38_legacy_variants')

    def test_normalize_reports_vcf(self):
        file_reports = aggregate_reports.normalize_reports(self.vcf_test_file, self.columns, self.genome_regions_symbol_dict)
        first_report = file_reports[0]
        self.assertEqual(len(file_reports), 2)
        self.assertEqual(first_report[aggregate_reports.COLUMN_SOURCE], '1000_Genomes')

    def test_normalize_reports_tsv(self):
        file_reports = aggregate_reports.normalize_reports(self.tsv_test_file, self.columns, self.genome_regions_symbol_dict)
        first_report = file_reports[0]
        self.assertEqual(len(file_reports), 2)
        self.assertEqual(first_report[aggregate_reports.COLUMN_SOURCE], 'ENIGMA')

    def test_get_reports_files(self):
        reports_files = aggregate_reports.get_reports_files(INPUT_DIRECTORY)
        self.assertEqual(len(reports_files), 10)
        self.assertNotIn("1000_Genomesready.vcf", reports_files)

    def test_aggregate_reports(self):
        reports_files = [INPUT_DIRECTORY + r for r in aggregate_reports.get_reports_files(INPUT_DIRECTORY)]
        reports = aggregate_reports.aggregate_reports(reports_files, self.columns, self.genome_regions_symbol_dict)

        # Each test file contains two reports, check that all reports are present
        self.assertEqual(len(reports), (len(self.sources) * 2))

        # Check that two of each source are present
        source_reports = {}
        source_column_index = self.columns.index("Source")
        for report in reports:
            source = report[source_column_index]
            if source not in source_reports:
                source_reports[source] = 1
            else:
                source_reports[source] += 1
        for source in self.sources:
            self.assertEqual(source_reports[source], 2)

    def test_aggregate_reports_maintains_proper_variant_effect_lovd_formatting(self):
        LOVD_reports_file = [INPUT_DIRECTORY + r for r in aggregate_reports.get_reports_files(INPUT_DIRECTORY) if r == 'LOVD.vcf']
        reports = aggregate_reports.aggregate_reports(LOVD_reports_file, self.columns, self.genome_regions_symbol_dict)

        # Check that two of each source are present
        variant_effect_lovd_index = self.columns.index("Variant_effect_LOVD")

        for variant in reports:
            self.assertIn(variant[variant_effect_lovd_index][0], ['?/.', '+/+'])



if __name__ == '__main__':
    pass
