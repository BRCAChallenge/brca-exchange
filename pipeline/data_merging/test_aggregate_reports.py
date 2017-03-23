import pytest
import unittest
import tempfile
import csv
from os import path
import aggregate_reports
import pdb


VCF_TESTDATA_FILENAME = path.join(path.dirname(__file__), '1000_Genomesready.vcf')
TSV_TESTDATA_FILENAME = path.join(path.dirname(__file__), 'ENIGMA_combined_with_bx_ids.tsv')


class TestStringMethods(unittest.TestCase):

    def setUp(self):
        self.columns = ['Source', 'Gene_symbol_ENIGMA', 'Genomic_Coordinate', 'Chr', 'Pos', 'Ref', 'Alt', 'Reference_sequence_ENIGMA', 'HGVS_cDNA_ENIGMA', 'BIC_Nomenclature_ENIGMA', 'Abbrev_AA_change_ENIGMA', 'URL_ENIGMA', 'Condition_ID_type_ENIGMA', 'Condition_ID_value_ENIGMA', 'Condition_category_ENIGMA', 'Clinical_significance_ENIGMA', 'Date_last_evaluated_ENIGMA', 'Assertion_method_ENIGMA', 'Assertion_method_citation_ENIGMA', 'Clinical_significance_citations_ENIGMA', 'Comment_on_clinical_significance_ENIGMA', 'Collection_method_ENIGMA', 'Allele_origin_ENIGMA', 'ClinVarAccession_ENIGMA', 'HGVS_protein_ENIGMA', 'BX_ID_ENIGMA', 'Clinical_Significance_ClinVar', 'Date_Last_Updated_ClinVar', 'BX_ID_ClinVar', 'HGVS_ClinVar', 'Submitter_ClinVar', 'Protein_ClinVar', 'SCV_ClinVar', 'Allele_Origin_ClinVar', 'Method_ClinVar', 'BX_ID_LOVD', 'Variant_frequency_LOVD', 'HGVS_cDNA_LOVD', 'HGVS_protein_LOVD', 'BX_ID_ESP', 'Minor_allele_frequency_ESP', 'polyPhen2_result_ESP', 'EUR_Allele_frequency_1000_Genomes', 'AFR_Allele_frequency_1000_Genomes', 'AMR_Allele_frequency_1000_Genomes', 'EAS_Allele_frequency_1000_Genomes', 'BX_ID_1000_Genomes', 'Allele_frequency_1000_Genomes', 'SAS_Allele_frequency_1000_Genomes', 'Allele_frequency_ExAC', 'BX_ID_ExAC', 'BX_ID_BIC', 'Patient_nationality_BIC', 'Clinical_importance_BIC', 'Clinical_classification_BIC', 'BIC_Designation_BIC', 'Literature_citation_BIC', 'Number_of_family_member_carrying_mutation_BIC', 'Germline_or_Somatic_BIC', 'Ethnicity_BIC', 'Mutation_type_BIC', 'IARC_class_exLOVD', 'BIC_Nomenclature_exLOVD', 'Sum_family_LR_exLOVD', 'Combined_prior_probablility_exLOVD', 'BX_ID_exLOVD', 'HGVS_cDNA_exLOVD', 'Literature_source_exLOVD', 'Co_occurrence_LR_exLOVD', 'Posterior_probability_exLOVD', 'Missense_analysis_prior_probability_exLOVD', 'Segregation_LR_exLOVD', 'HGVS_protein_exLOVD']
        self.vcf_test_file = VCF_TESTDATA_FILENAME
        self.tsv_test_file = TSV_TESTDATA_FILENAME

    def test_normalize_reports_vcf(self):
        file_reports = aggregate_reports.normalize_reports(self.vcf_test_file, self.columns)
        first_report = file_reports[0]
        self.assertEqual(len(file_reports), 3)
        self.assertEqual(first_report[aggregate_reports.COLUMN_SOURCE], '1000_Genomes')

    def test_normalize_reports_tsv(self):
        file_reports = aggregate_reports.normalize_reports(self.tsv_test_file, self.columns)
        first_report = file_reports[0]
        self.assertEqual(len(file_reports), 2)
        self.assertEqual(first_report[aggregate_reports.COLUMN_SOURCE], 'ENIGMA')


if __name__ == '__main__':
    pass
