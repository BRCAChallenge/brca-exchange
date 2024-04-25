import pytest
import unittest
from .aggregate_across_columns import selectMaxAlleleFrequency, selectAlleleFrequency, FIELDS_TO_REMOVE, FIELDS_TO_ADD, FIELDS_TO_RENAME, setOutputColumns, update_basic_fields, EMPTY

class TestStringMethods(unittest.TestCase):

    def setUp(self):

        self.oldRow = {
            'Variant_frequency_LOVD': EMPTY,
            'ClinVarAccession_ENIGMA': 'SCV000244909',
            'Chr': '13',
            'Abbrev_AA_change_ENIGMA': EMPTY,
            'BIC_Nomenclature_exLOVD': EMPTY,
            'Co_occurrence_LR_exLOVD': EMPTY,
            'Submitter_ClinVar': 'Evidence-based_Network_for_the_Interpretation_of_Germline_Mutant_Alleles_(ENIGMA)',
            'Submitters_LOVD': EMPTY,
            'Method_ClinVar': 'curation',
            'Pathogenicity_all': 'Benign(ENIGMA); Benign (ClinVar)',
            'BIC_Nomenclature_ENIGMA': '-536A>G',
            'Literature_source_exLOVD': EMPTY,
            'Collection_method_ENIGMA': 'Curation',
            'Sum_family_LR_exLOVD': EMPTY,
            'HGVS_cDNA_LOVD': EMPTY,
            'Individuals_LOVD': EMPTY,
            'URL_ENIGMA': EMPTY,
            'BX_ID_exLOVD': EMPTY,
            'Allele_Origin_ClinVar': 'germline',
            'Source': 'ENIGMA,ClinVar,1000_Genomes',
            'Condition_ID_value_ENIGMA': 'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)',
            'HGVS_protein_ENIGMA': 'p.?',
            'Ref': 'A',
            'BX_ID_LOVD': EMPTY,
            'Gene_symbol_ENIGMA': EMPTY,
            'Comment_on_clinical_significance_ENIGMA': 'Class 1 not pathogenic based on frequency >1% in an outbred sampleset. Frequency 0.14 (African), derived from 1000 genomes (2012-04-30).',
            'Missense_analysis_prior_probability_exLOVD': EMPTY,
            'Assertion_method_ENIGMA': 'ENIGMA BRCA1/2 Classification Criteria (2015)',
            'Posterior_probability_exLOVD': EMPTY,
            'Polyphen_Score': EMPTY,
            'Reference_sequence_ENIGMA': 'NM_000059.3',
            'HGVS_cDNA_exLOVD': EMPTY,
            'Hg38_End': EMPTY,
            'HGVS_cDNA_ENIGMA': 'c.-764A>G',
            'RNA_LOVD': EMPTY,
            'BX_ID_ClinVar': '1',
            'IARC_class_exLOVD': EMPTY,
            'Allele_origin_ENIGMA': 'Germline',
            'Date_last_evaluated_ENIGMA': '1/12/2015',
            'HGVS_ClinVar': 'NM_000059.3.c.-764A>G',
            'Date_Last_Updated_ClinVar': '2015-01-12',
            'SCV_ClinVar': 'SCV000244909',
            'Pathogenicity_expert': 'Benign / Little Clinical Significance',
            'Combined_prior_probablility_exLOVD': EMPTY,
            'Condition_category_ENIGMA': 'Disease',
            'Genetic_origin_LOVD': EMPTY,
            'Segregation_LR_exLOVD': EMPTY,
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Assertion_method_citation_ENIGMA': 'https://enigmaconsortium.org/library/general-documents/',
            'Condition_ID_type_ENIGMA': 'OMIM',
            'HGVS_protein_LOVD': EMPTY,
            'Clinical_significance_ENIGMA': 'Benign',
            'Clinical_Significance_ClinVar': 'Benign',
            'BX_ID_ENIGMA': '2890',
            'Sift_Prediction': EMPTY,
            'HGVS_RNA': EMPTY,
            'Clinical_significance_citations_ENIGMA': EMPTY,
            'Variant_effect_LOVD': EMPTY,
            'Polyphen_Prediction': EMPTY,
            'Sift_Score': EMPTY,
            'Genomic_Coordinate': 'chr13:32314943:A>G',
            'Alt': 'G',
            'Protein_ClinVar': 'None',
            'Hg38_Start': EMPTY,
            'Pos': '32314943',
            'HGVS_protein_exLOVD': EMPTY
        }

        self.newRow = {
            'Variant_frequency_LOVD': EMPTY,
            'ClinVarAccession_ENIGMA': 'SCV000244909',
            'Chr': '13',
            'Protein_Change': EMPTY,
            'BIC_Nomenclature_exLOVD': EMPTY,
            'Co_occurrence_LR_exLOVD': EMPTY,
            'Submitter_ClinVar': 'Evidence-based_Network_for_the_Interpretation_of_Germline_Mutant_Alleles_(ENIGMA)',
            'Submitters_LOVD': EMPTY,
            'Method_ClinVar': 'curation',
            'Pathogenicity_all': 'Benign(ENIGMA); Benign (ClinVar)',
            'BIC_Nomenclature': '-536A>G',
            'Literature_source_exLOVD': EMPTY,
            'Collection_method_ENIGMA': 'Curation',
            'Sum_family_LR_exLOVD': EMPTY,
            'HGVS_cDNA_LOVD': EMPTY,
            'Individuals_LOVD': EMPTY,
            'URL_ENIGMA': EMPTY,
            'BX_ID_exLOVD': EMPTY,
            'Allele_Origin_ClinVar': 'germline',
            'Source': 'ENIGMA,ClinVar',
            'Condition_ID_value_ENIGMA': 'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)',
            'HGVS_Protein': 'p.?',
            'Ref': 'A',
            'BX_ID_LOVD': EMPTY,
            'Gene_Symbol': 'BRCA2',
            'Comment_on_clinical_significance_ENIGMA': 'Class 1 not pathogenic based on frequency >1% in an outbred sampleset. Frequency 0.14 (African), derived from 1000 genomes (2012-04-30).',
            'Missense_analysis_prior_probability_exLOVD': EMPTY,
            'Assertion_method_ENIGMA': 'ENIGMA BRCA1/2 Classification Criteria (2015)',
            'Posterior_probability_exLOVD': EMPTY,
            'Polyphen_Score': EMPTY,
            'Reference_Sequence': 'NM_000059.3',
            'HGVS_cDNA_exLOVD': EMPTY,
            'Hg38_End': 32314943, 'HGVS_cDNA': 'c.-764A>G',
            'SAS_Allele_frequency_1000_Genomes': '0.0',
            'RNA_LOVD': EMPTY,
            'BX_ID_ClinVar': '1',
            'IARC_class_exLOVD': EMPTY,
            'Allele_origin_ENIGMA': 'Germline',
            'Date_last_evaluated_ENIGMA': '1/12/2015',
            'HGVS_ClinVar': 'NM_000059.3.c.-764A>G',
            'Date_Last_Updated_ClinVar': '2015-01-12',
            'SCV_ClinVar': 'SCV000244909',
            'Pathogenicity_expert': 'Benign / Little Clinical Significance',
            'Allele_frequency_1000_Genomes': '0.0401358',
            'Combined_prior_probablility_exLOVD': EMPTY,
            'AMR_Allele_frequency_1000_Genomes': '0.0072',
            'Condition_category_ENIGMA': 'Disease',
            'Genetic_origin_LOVD': EMPTY,
            'Segregation_LR_exLOVD': EMPTY,
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Assertion_method_citation_ENIGMA': 'https://enigmaconsortium.org/library/general-documents/',
            'Condition_ID_type_ENIGMA': 'OMIM',
            'HGVS_protein_LOVD': EMPTY,
            'Clinical_significance_ENIGMA': 'Benign',
            'Clinical_Significance_ClinVar': 'Benign',
            'BX_ID_ENIGMA': '2890',
            'Sift_Prediction': EMPTY,
            'Clinical_significance_citations_ENIGMA': EMPTY,
            'Variant_effect_LOVD': EMPTY,
            'Polyphen_Prediction': EMPTY,
            'Sift_Score': EMPTY,
            'Genomic_Coordinate_hg38': 'chr13:32314943:A>G',
            'Alt': 'G',
            'Protein_ClinVar': 'None',
            'Hg38_Start': '32314943',
            'Pos': '32314943',
            'HGVS_protein_exLOVD': EMPTY
        }

        self.initialFields = list(FIELDS_TO_RENAME.keys()) + FIELDS_TO_REMOVE

        self.newRowAlleleFrequencies = {
            'Allele_frequency_genome_GnomAD': EMPTY,
            'Allele_frequency_exome_GnomAD': EMPTY,
            'Allele_count_genome_GnomAD': EMPTY,
            'Allele_count_exome_GnomAD': EMPTY,
            'Allele_number_genome_GnomAD': EMPTY,
            'Allele_number_exome_GnomAD': EMPTY,
        }


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

        for oldName, newName in FIELDS_TO_RENAME.items():
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

        for oldName, newName in FIELDS_TO_RENAME.items():
            self.assertNotIn(oldName, list(updatedRow.keys()))
            self.assertIn(newName, list(updatedRow.keys()))

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
        for attr, value in self.newRowAlleleFrequencies.items():
            self.newRowAlleleFrequencies[attr] = '-'


    def test_determine_gnomAD_allele_frequency(self):
        self.newRowAlleleFrequencies['Allele_frequency_genome_GnomAD'] = '0.345'
        self.newRowAlleleFrequencies['Allele_count_genome_GnomAD'] = '345'
        self.newRowAlleleFrequencies['Allele_count_exome_GnomAD'] = '0'
        self.newRowAlleleFrequencies['Allele_number_genome_GnomAD'] = '1000'
        self.newRowAlleleFrequencies['Allele_number_exome_GnomAD'] = '0'

        AF = selectAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEqual(AF, '0.345 (GnomAD)')

        self.newRowAlleleFrequencies['Allele_frequency_exome_GnomAD'] = '0.5'
        self.newRowAlleleFrequencies['Allele_count_genome_GnomAD'] = '1'
        self.newRowAlleleFrequencies['Allele_count_exome_GnomAD'] = '1'
        self.newRowAlleleFrequencies['Allele_number_genome_GnomAD'] = '2'
        self.newRowAlleleFrequencies['Allele_number_exome_GnomAD'] = '2'

        AF = selectAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEqual(AF, '0.5 (GnomAD)')


if __name__ == '__main__':
    pass
