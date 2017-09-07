import pytest
import unittest
from aggregate_across_columns import selectMaxAlleleFrequency, selectAlleleFrequency, FIELDS_TO_REMOVE, FIELDS_TO_ADD, FIELDS_TO_RENAME, setOutputColumns, update_basic_fields, EMPTY

class TestStringMethods(unittest.TestCase):

    def setUp(self):

        self.oldRow = {
            'Variant_frequency_LOVD': EMPTY,
            'Allele_frequency_FIN_ExAC': EMPTY,
            'ClinVarAccession_ENIGMA': 'SCV000244909',
            'Homozygous_count_AFR_ExAC': EMPTY,
            'BX_ID_ExAC': EMPTY,
            'Allele_frequency_AFR_ExAC': EMPTY,
            'Chr': '13',
            'Abbrev_AA_change_ENIGMA': EMPTY,
            'BIC_Nomenclature_exLOVD': EMPTY,
            'Co_occurrence_LR_exLOVD': EMPTY,
            'Homozygous_count_EAS_ExAC': EMPTY,
            'Submitter_ClinVar': 'Evidence-based_Network_for_the_Interpretation_of_Germline_Mutant_Alleles_(ENIGMA)',
            'Allele_frequency_EAS_ExAC': EMPTY,
            'Submitters_LOVD': EMPTY,
            'Clinical_classification_BIC': EMPTY,
            'Homozygous_count_NFE_ExAC': EMPTY,
            'Allele_count_SAS_ExAC': EMPTY,
            'Method_ClinVar': 'curation',
            'Allele_count_NFE_ExAC': EMPTY,
            'Pathogenicity_all': 'Benign(ENIGMA); Benign (ClinVar)',
            'Germline_or_Somatic_BIC': EMPTY,
            'Homozygous_count_SAS_ExAC': EMPTY,
            'BIC_Nomenclature_ENIGMA': '-536A>G',
            'Allele_number_FIN_ExAC': EMPTY,
            'Literature_source_exLOVD': EMPTY,
            'Collection_method_ENIGMA': 'Curation',
            'Sum_family_LR_exLOVD': EMPTY,
            'HGVS_cDNA_LOVD': EMPTY,
            'Homozygous_count_FIN_ExAC': EMPTY,
            'EAS_Allele_frequency_1000_Genomes': '0.0',
            'Ethnicity_BIC': EMPTY,
            'Individuals_LOVD': EMPTY,
            'polyPhen2_result_ESP': EMPTY,
            'URL_ENIGMA': EMPTY,
            'BX_ID_exLOVD': EMPTY,
            'Allele_Origin_ClinVar': 'germline',
            'Allele_frequency_AMR_ExAC': EMPTY,
            'AFR_Allele_frequency_1000_Genomes': '0.1475',
            'EUR_Allele_frequency_1000_Genomes': '0.001',
            'Source': 'ENIGMA,ClinVar,1000_Genomes',
            'Condition_ID_value_ENIGMA': 'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)',
            'HGVS_protein_ENIGMA': 'p.?',
            'Ref': 'A',
            'Allele_number_AFR_ExAC': EMPTY,
            'Allele_count_AFR_ExAC': EMPTY,
            'BX_ID_LOVD': EMPTY,
            'Gene_symbol_ENIGMA': EMPTY,
            'Comment_on_clinical_significance_ENIGMA': 'Class 1 not pathogenic based on frequency >1% in an outbred sampleset. Frequency 0.14 (African), derived from 1000 genomes (2012-04-30).',
            'Missense_analysis_prior_probability_exLOVD': EMPTY,
            'Assertion_method_ENIGMA': 'ENIGMA BRCA1/2 Classification Criteria (2015)',
            'Posterior_probability_exLOVD': EMPTY,
            'Polyphen_Score': EMPTY,
            'Reference_sequence_ENIGMA': 'NM_000059.3',
            'HGVS_cDNA_exLOVD': EMPTY,
            'Allele_count_EAS_ExAC': EMPTY,
            'Hg38_End': EMPTY,
            'HGVS_cDNA_ENIGMA': 'c.-764A>G',
            'SAS_Allele_frequency_1000_Genomes': '0.0',
            'RNA_LOVD': EMPTY,
            'BX_ID_1000_Genomes': '1',
            'BX_ID_ClinVar': '1',
            'IARC_class_exLOVD': EMPTY,
            'BX_ID_BIC': EMPTY,
            'Allele_number_NFE_ExAC': EMPTY,
            'Allele_origin_ENIGMA': 'Germline',
            'Allele_number_OTH_ExAC': EMPTY,
            'Date_last_evaluated_ENIGMA': '1/12/2015',
            'HGVS_ClinVar': 'NM_000059.3.c.-764A>G',
            'BIC_Designation_BIC': EMPTY,
            'Allele_frequency_SAS_ExAC': EMPTY,
            'Date_Last_Updated_ClinVar': '2015-01-12',
            'Allele_number_EAS_ExAC': EMPTY,
            'Allele_frequency_OTH_ExAC': EMPTY,
            'SCV_ClinVar': 'SCV000244909',
            'Pathogenicity_expert': 'Benign / Little Clinical Significance',
            'Allele_frequency_1000_Genomes': '0.0401358',
            'Combined_prior_probablility_exLOVD': EMPTY,
            'AMR_Allele_frequency_1000_Genomes': '0.0072',
            'Condition_category_ENIGMA': 'Disease',
            'Allele_count_AMR_ExAC': EMPTY,
            'BX_ID_ESP': EMPTY,
            'Patient_nationality_BIC': EMPTY,
            'Genetic_origin_LOVD': EMPTY,
            'Number_of_family_member_carrying_mutation_BIC': EMPTY,
            'Segregation_LR_exLOVD': EMPTY,
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Allele_frequency_ExAC': EMPTY,
            'Mutation_type_BIC': EMPTY,
            'Assertion_method_citation_ENIGMA': 'https://enigmaconsortium.org/library/general-documents/',
            'Condition_ID_type_ENIGMA': 'OMIM',
            'Allele_count_OTH_ExAC': EMPTY,
            'HGVS_protein_LOVD': EMPTY,
            'Clinical_importance_BIC': EMPTY,
            'Homozygous_count_OTH_ExAC': EMPTY,
            'Allele_count_FIN_ExAC': EMPTY,
            'Clinical_significance_ENIGMA': 'Benign',
            'Minor_allele_frequency_percent_ESP': EMPTY,
            'Allele_Frequency_ESP': EMPTY,
            'Homozygous_count_AMR_ExAC': EMPTY,
            'Clinical_Significance_ClinVar': 'Benign',
            'AA_Allele_Frequency_ESP': EMPTY,
            'BX_ID_ENIGMA': '2890',
            'Sift_Prediction': EMPTY,
            'EA_Allele_Frequency_ESP': EMPTY,
            'HGVS_RNA': EMPTY,
            'Clinical_significance_citations_ENIGMA': EMPTY,
            'Variant_effect_LOVD': EMPTY,
            'Polyphen_Prediction': EMPTY,
            'Sift_Score': EMPTY,
            'Genomic_Coordinate': 'chr13:32314943:A>G',
            'Alt': 'G',
            'Protein_ClinVar': 'None',
            'Literature_citation_BIC': EMPTY,
            'Allele_frequency_NFE_ExAC': EMPTY,
            'Hg38_Start': EMPTY,
            'Pos': '32314943',
            'HGVS_protein_exLOVD': EMPTY,
            'Allele_number_SAS_ExAC': EMPTY,
            'Allele_number_AMR_ExAC': EMPTY
        }

        self.newRow = {
            'Variant_frequency_LOVD': EMPTY,
            'Allele_frequency_FIN_ExAC': EMPTY,
            'ClinVarAccession_ENIGMA': 'SCV000244909',
            'Homozygous_count_AFR_ExAC': EMPTY,
            'BX_ID_ExAC': EMPTY,
            'Allele_frequency_AFR_ExAC': EMPTY,
            'Chr': '13',
            'Protein_Change': EMPTY,
            'BIC_Nomenclature_exLOVD': EMPTY,
            'Co_occurrence_LR_exLOVD': EMPTY,
            'Homozygous_count_EAS_ExAC': EMPTY,
            'Submitter_ClinVar': 'Evidence-based_Network_for_the_Interpretation_of_Germline_Mutant_Alleles_(ENIGMA)',
            'Allele_frequency_EAS_ExAC': EMPTY,
            'Submitters_LOVD': EMPTY,
            'Clinical_classification_BIC': EMPTY,
            'Homozygous_count_NFE_ExAC': EMPTY,
            'Allele_count_SAS_ExAC': EMPTY,
            'Method_ClinVar': 'curation',
            'Allele_count_NFE_ExAC': EMPTY,
            'Pathogenicity_all': 'Benign(ENIGMA); Benign (ClinVar)',
            'Germline_or_Somatic_BIC': EMPTY,
            'Homozygous_count_SAS_ExAC': EMPTY,
            'BIC_Nomenclature': '-536A>G',
            'Allele_number_FIN_ExAC': EMPTY,
            'Literature_source_exLOVD': EMPTY,
            'Collection_method_ENIGMA': 'Curation',
            'Sum_family_LR_exLOVD': EMPTY,
            'HGVS_cDNA_LOVD': EMPTY,
            'Homozygous_count_FIN_ExAC': EMPTY,
            'EAS_Allele_frequency_1000_Genomes': '0.0',
            'Ethnicity_BIC': EMPTY,
            'Individuals_LOVD': EMPTY,
            'polyPhen2_result_ESP': EMPTY,
            'URL_ENIGMA': EMPTY,
            'BX_ID_exLOVD': EMPTY,
            'Allele_Origin_ClinVar': 'germline',
            'Allele_frequency_AMR_ExAC': EMPTY,
            'AFR_Allele_frequency_1000_Genomes': '0.1475',
            'EUR_Allele_frequency_1000_Genomes': '0.001',
            'Source': 'ENIGMA,ClinVar,1000_Genomes',
            'Condition_ID_value_ENIGMA': 'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)',
            'HGVS_Protein': 'p.?',
            'Ref': 'A',
            'Allele_number_AFR_ExAC': EMPTY,
            'Allele_count_AFR_ExAC': EMPTY,
            'BX_ID_LOVD': EMPTY,
            'Gene_Symbol': 'BRCA2',
            'Comment_on_clinical_significance_ENIGMA': 'Class 1 not pathogenic based on frequency >1% in an outbred sampleset. Frequency 0.14 (African), derived from 1000 genomes (2012-04-30).',
            'Missense_analysis_prior_probability_exLOVD': EMPTY,
            'Assertion_method_ENIGMA': 'ENIGMA BRCA1/2 Classification Criteria (2015)',
            'Posterior_probability_exLOVD': EMPTY,
            'Polyphen_Score': EMPTY,
            'Reference_Sequence': 'NM_000059.3',
            'HGVS_cDNA_exLOVD': EMPTY,
            'Allele_count_EAS_ExAC': EMPTY,
            'Hg38_End': 32314943, 'HGVS_cDNA': 'c.-764A>G',
            'SAS_Allele_frequency_1000_Genomes': '0.0',
            'RNA_LOVD': EMPTY,
            'BX_ID_1000_Genomes': '1',
            'BX_ID_ClinVar': '1',
            'IARC_class_exLOVD': EMPTY,
            'BX_ID_BIC': EMPTY,
            'Allele_number_NFE_ExAC': EMPTY,
            'Allele_origin_ENIGMA': 'Germline',
            'Allele_number_OTH_ExAC': EMPTY,
            'Date_last_evaluated_ENIGMA': '1/12/2015',
            'HGVS_ClinVar': 'NM_000059.3.c.-764A>G',
            'BIC_Designation_BIC': EMPTY,
            'Allele_frequency_SAS_ExAC': EMPTY,
            'Date_Last_Updated_ClinVar': '2015-01-12',
            'Allele_number_EAS_ExAC': EMPTY,
            'Allele_frequency_OTH_ExAC': EMPTY,
            'SCV_ClinVar': 'SCV000244909',
            'Pathogenicity_expert': 'Benign / Little Clinical Significance',
            'Allele_frequency_1000_Genomes': '0.0401358',
            'Combined_prior_probablility_exLOVD': EMPTY,
            'AMR_Allele_frequency_1000_Genomes': '0.0072',
            'Condition_category_ENIGMA': 'Disease',
            'Allele_count_AMR_ExAC': EMPTY,
            'BX_ID_ESP': EMPTY,
            'Patient_nationality_BIC': EMPTY,
            'Genetic_origin_LOVD': EMPTY,
            'Number_of_family_member_carrying_mutation_BIC': EMPTY,
            'Segregation_LR_exLOVD': EMPTY,
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Allele_frequency_ExAC': EMPTY,
            'Mutation_type_BIC': EMPTY,
            'Assertion_method_citation_ENIGMA': 'https://enigmaconsortium.org/library/general-documents/',
            'Condition_ID_type_ENIGMA': 'OMIM',
            'Allele_count_OTH_ExAC': EMPTY,
            'HGVS_protein_LOVD': EMPTY,
            'Clinical_importance_BIC': EMPTY,
            'Homozygous_count_OTH_ExAC': EMPTY,
            'Allele_count_FIN_ExAC': EMPTY,
            'Clinical_significance_ENIGMA': 'Benign',
            'Minor_allele_frequency_percent_ESP': EMPTY,
            'Allele_Frequency_ESP': EMPTY,
            'Homozygous_count_AMR_ExAC': EMPTY,
            'Clinical_Significance_ClinVar': 'Benign',
            'AA_Allele_Frequency_ESP': EMPTY,
            'BX_ID_ENIGMA': '2890',
            'Sift_Prediction': EMPTY,
            'EA_Allele_Frequency_ESP': EMPTY,
            'Clinical_significance_citations_ENIGMA': EMPTY,
            'Variant_effect_LOVD': EMPTY,
            'Polyphen_Prediction': EMPTY,
            'Sift_Score': EMPTY,
            'Genomic_Coordinate_hg38': 'chr13:32314943:A>G',
            'Alt': 'G',
            'Protein_ClinVar': 'None',
            'Literature_citation_BIC': EMPTY,
            'Allele_frequency_NFE_ExAC': EMPTY,
            'Hg38_Start': '32314943',
            'Pos': '32314943',
            'HGVS_protein_exLOVD': EMPTY,
            'Allele_number_SAS_ExAC': EMPTY,
            'Allele_number_AMR_ExAC': EMPTY
        }

        self.initialFields = FIELDS_TO_RENAME.keys() + FIELDS_TO_REMOVE

        self.newRowAlleleFrequencies = {
            'Allele_frequency_FIN_ExAC': EMPTY,
            'Allele_frequency_AFR_ExAC': EMPTY,
            'Allele_frequency_EAS_ExAC': EMPTY,
            'EAS_Allele_frequency_1000_Genomes': '0.0',
            'Allele_frequency_AMR_ExAC': EMPTY,
            'AFR_Allele_frequency_1000_Genomes': '0.1475',
            'EUR_Allele_frequency_1000_Genomes': '0.001',
            'SAS_Allele_frequency_1000_Genomes': '0.0',
            'Allele_frequency_SAS_ExAC': EMPTY,
            'Allele_frequency_OTH_ExAC': EMPTY,
            'Allele_frequency_1000_Genomes': '0.0401358',
            'AMR_Allele_frequency_1000_Genomes': '0.0072',
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Allele_frequency_ExAC': EMPTY,
            'Minor_allele_frequency_percent_ESP': EMPTY,
            'Allele_Frequency_ESP': EMPTY,
            'AA_Allele_Frequency_ESP': EMPTY,
            'EA_Allele_Frequency_ESP': EMPTY,
            'Allele_frequency_NFE_ExAC': EMPTY,
        }

    def test_select_max_allele_frequency(self):
        maxFreqString = selectMaxAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(maxFreqString, EMPTY)

        self.newRowAlleleFrequencies["Allele_frequency_SAS_ExAC"] = '0.305'
        maxFreqString = selectMaxAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(maxFreqString, EMPTY)

        for attr, value in self.newRowAlleleFrequencies.iteritems():
            self.newRowAlleleFrequencies[attr] = EMPTY

        maxFreqString = selectMaxAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(maxFreqString, EMPTY)

    def test_set_output_columns(self):
        '''
        Tests that: 
        none of the fields in FIELDS_TO_REMOVE are present, 
        all of the fields from FIELDS_TO_ADD are present, 
        and that all of the fields in FIEDLS_TO_RENAME have been properly renamed.
        '''
        outputFields = setOutputColumns(self.initialFields, FIELDS_TO_REMOVE, FIELDS_TO_ADD, FIELDS_TO_RENAME)

        for field in FIELDS_TO_REMOVE:
            self.assertNotIn(field, outputFields)
            
        for field in FIELDS_TO_ADD:
            self.assertIn(field, outputFields)
            
        for oldName, newName in FIELDS_TO_RENAME.iteritems():
            self.assertNotIn(oldName, outputFields)
            self.assertIn(newName, outputFields)

    def test_update_basic_fields(self):
        '''
        Tests that:
        all fields have been appropriately renamed,
        Hg38_Start and Hg38_End are set to the correct values,
        Gene_Symbol is set correctly based on the chromosome number,
        and HGVS_RNA is set to the correct value.
        '''
        oldRow_copy1 = self.oldRow.copy()
        oldRow_copy2 = self.oldRow.copy()
        updatedRow = update_basic_fields(self.oldRow, FIELDS_TO_RENAME)
        
        for oldName, newName in FIELDS_TO_RENAME.iteritems():
            self.assertNotIn(oldName, updatedRow.keys())
            self.assertIn(newName, updatedRow.keys())

        expected_start = self.oldRow["Pos"]
        self.assertEqual(expected_start, updatedRow["Hg38_Start"])

        expected_end = int(expected_start) + len(self.oldRow["Ref"]) - 1
        self.assertEqual(expected_end, updatedRow["Hg38_End"])

        oldRow_copy1["Genomic_Coordinate"] = 'chr17:32314943:A>G'
        updatedRow = update_basic_fields(oldRow_copy1, FIELDS_TO_RENAME)
        self.assertEqual(updatedRow["Gene_Symbol"], 'BRCA1')

        oldRow_copy2["Genomic_Coordinate"] = 'chr13:32314943:A>G'
        updatedRow = update_basic_fields(oldRow_copy2, FIELDS_TO_RENAME)
        self.assertEqual(updatedRow["Gene_Symbol"], 'BRCA2')

        self.assertEqual(updatedRow["HGVS_RNA"], EMPTY)

    def test_select_allele_frequency(self):
        for attr, value in self.newRowAlleleFrequencies.iteritems():
            self.newRowAlleleFrequencies[attr] = '-'


        self.newRowAlleleFrequencies['Minor_allele_frequency_percent_ESP'] = '20.345'
        AF = selectAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(AF, '0.20345 (ESP)')

        self.newRowAlleleFrequencies['Minor_allele_frequency_percent_ESP'] = '0.0'
        AF = selectAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(AF, '0.0 (ESP)')

        self.newRowAlleleFrequencies['Minor_allele_frequency_percent_ESP'] = '2'
        AF = selectAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(AF, '0.02 (ESP)')


if __name__ == '__main__':
    pass
