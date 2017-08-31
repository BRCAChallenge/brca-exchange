import pytest
import unittest
from aggregate_across_columns import selectMaxAlleleFrequency, selectAlleleFrequency


class TestStringMethods(unittest.TestCase):

    def setUp(self):
        self.newRow = {
            'Variant_frequency_LOVD': '-',
            'Allele_frequency_FIN_ExAC': '-',
            'ClinVarAccession_ENIGMA': 'SCV000244909',
            'Homozygous_count_AFR_ExAC': '-',
            'BX_ID_ExAC': '-',
            'Allele_frequency_AFR_ExAC': '-',
            'Chr': '13',
            'Protein_Change': '-',
            'BIC_Nomenclature_exLOVD': '-',
            'Co_occurrence_LR_exLOVD': '-',
            'Homozygous_count_EAS_ExAC': '-',
            'Submitter_ClinVar': 'Evidence-based_Network_for_the_Interpretation_of_Germline_Mutant_Alleles_(ENIGMA)',
            'Allele_frequency_EAS_ExAC': '-',
            'Submitters_LOVD': '-',
            'Clinical_classification_BIC': '-',
            'Homozygous_count_NFE_ExAC': '-',
            'Allele_count_SAS_ExAC': '-',
            'Method_ClinVar': 'curation',
            'Allele_count_NFE_ExAC': '-',
            'Pathogenicity_all': 'Benign(ENIGMA); Benign (ClinVar)',
            'Germline_or_Somatic_BIC': '-',
            'Homozygous_count_SAS_ExAC': '-',
            'BIC_Nomenclature': '-536A>G',
            'Allele_number_FIN_ExAC': '-',
            'Literature_source_exLOVD': '-',
            'Collection_method_ENIGMA': 'Curation',
            'Sum_family_LR_exLOVD': '-',
            'HGVS_cDNA_LOVD': '-',
            'Homozygous_count_FIN_ExAC': '-',
            'EAS_Allele_frequency_1000_Genomes': '0.0',
            'Ethnicity_BIC': '-',
            'Individuals_LOVD': '-',
            'polyPhen2_result_ESP': '-',
            'URL_ENIGMA': '-',
            'BX_ID_exLOVD': '-',
            'Allele_Origin_ClinVar': 'germline',
            'Allele_frequency_AMR_ExAC': '-',
            'AFR_Allele_frequency_1000_Genomes': '0.1475',
            'EUR_Allele_frequency_1000_Genomes': '0.001',
            'Source': 'ENIGMA,ClinVar,1000_Genomes',
            'Condition_ID_value_ENIGMA': 'BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (612555)',
            'HGVS_Protein': 'p.?',
            'Ref': 'A',
            'Allele_number_AFR_ExAC': '-',
            'Allele_count_AFR_ExAC': '-',
            'BX_ID_LOVD': '-',
            'Gene_Symbol': 'BRCA2',
            'Comment_on_clinical_significance_ENIGMA': 'Class 1 not pathogenic based on frequency >1% in an outbred sampleset. Frequency 0.14 (African), derived from 1000 genomes (2012-04-30).',
            'Missense_analysis_prior_probability_exLOVD': '-',
            'Assertion_method_ENIGMA': 'ENIGMA BRCA1/2 Classification Criteria (2015)',
            'Posterior_probability_exLOVD': '-',
            'Polyphen_Score': '-',
            'Reference_Sequence': 'NM_000059.3',
            'HGVS_cDNA_exLOVD': '-',
            'Allele_count_EAS_ExAC': '-',
            'Hg38_End': 32314943, 'HGVS_cDNA': 'c.-764A>G',
            'SAS_Allele_frequency_1000_Genomes': '0.0',
            'RNA_LOVD': '-',
            'BX_ID_1000_Genomes': '1',
            'BX_ID_ClinVar': '1',
            'IARC_class_exLOVD': '-',
            'BX_ID_BIC': '-',
            'Allele_number_NFE_ExAC': '-',
            'Allele_origin_ENIGMA': 'Germline',
            'Allele_number_OTH_ExAC': '-',
            'Date_last_evaluated_ENIGMA': '1/12/2015',
            'HGVS_ClinVar': 'NM_000059.3.c.-764A>G',
            'BIC_Designation_BIC': '-',
            'Allele_frequency_SAS_ExAC': '-',
            'Date_Last_Updated_ClinVar': '2015-01-12',
            'Allele_number_EAS_ExAC': '-',
            'Allele_frequency_OTH_ExAC': '-',
            'SCV_ClinVar': 'SCV000244909',
            'Pathogenicity_expert': 'Benign / Little Clinical Significance',
            'Allele_frequency_1000_Genomes': '0.0401358',
            'Combined_prior_probablility_exLOVD': '-',
            'AMR_Allele_frequency_1000_Genomes': '0.0072',
            'Condition_category_ENIGMA': 'Disease',
            'Allele_count_AMR_ExAC': '-',
            'BX_ID_ESP': '-',
            'Patient_nationality_BIC': '-',
            'Genetic_origin_LOVD': '-',
            'Number_of_family_member_carrying_mutation_BIC': '-',
            'Segregation_LR_exLOVD': '-',
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Allele_frequency_ExAC': '-',
            'Mutation_type_BIC': '-',
            'Assertion_method_citation_ENIGMA': 'https://enigmaconsortium.org/library/general-documents/',
            'Condition_ID_type_ENIGMA': 'OMIM',
            'Allele_count_OTH_ExAC': '-',
            'HGVS_protein_LOVD': '-',
            'Clinical_importance_BIC': '-',
            'Homozygous_count_OTH_ExAC': '-',
            'Allele_count_FIN_ExAC': '-',
            'Clinical_significance_ENIGMA': 'Benign',
            'Minor_allele_frequency_percent_ESP': '-',
            'Allele_Frequency_ESP': '-',
            'Homozygous_count_AMR_ExAC': '-',
            'Clinical_Significance_ClinVar': 'Benign',
            'AA_Allele_Frequency_ESP': '-',
            'BX_ID_ENIGMA': '2890',
            'Sift_Prediction': '-',
            'EA_Allele_Frequency_ESP': '-',
            'HGVS_RNA': '-',
            'Clinical_significance_citations_ENIGMA': '-',
            'Variant_effect_LOVD': '-',
            'Polyphen_Prediction': '-',
            'Sift_Score': '-',
            'Genomic_Coordinate_hg38': 'chr13:32314943:A>G',
            'Alt': 'G',
            'Protein_ClinVar': 'None',
            'Literature_citation_BIC': '-',
            'Allele_frequency_NFE_ExAC': '-',
            'Hg38_Start': '32314943',
            'Pos': '32314943',
            'HGVS_protein_exLOVD': '-',
            'Allele_number_SAS_ExAC': '-',
            'Allele_number_AMR_ExAC': '-'
        }

        self.newRowAlleleFrequencies = {
            'Allele_frequency_FIN_ExAC': '-',
            'Allele_frequency_AFR_ExAC': '-',
            'Allele_frequency_EAS_ExAC': '-',
            'EAS_Allele_frequency_1000_Genomes': '0.0',
            'Allele_frequency_AMR_ExAC': '-',
            'AFR_Allele_frequency_1000_Genomes': '0.1475',
            'EUR_Allele_frequency_1000_Genomes': '0.001',
            'SAS_Allele_frequency_1000_Genomes': '0.0',
            'Allele_frequency_SAS_ExAC': '-',
            'Allele_frequency_OTH_ExAC': '-',
            'Allele_frequency_1000_Genomes': '0.0401358',
            'AMR_Allele_frequency_1000_Genomes': '0.0072',
            'Allele_Frequency': '0.0401358 (1000 Genomes)',
            'Allele_frequency_ExAC': '-',
            'Minor_allele_frequency_percent_ESP': '20.345',
            'Allele_Frequency_ESP': '-',
            'AA_Allele_Frequency_ESP': '-',
            'EA_Allele_Frequency_ESP': '-',
            'Allele_frequency_NFE_ExAC': '-',
        }

        self.maxAlleleFrequencyField = "AFR_Allele_frequency_1000_Genomes"

    def test_select_max_allele_frequency(self):
        maxFreqString = selectMaxAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(maxFreqString, '-')

        self.newRowAlleleFrequencies["Allele_frequency_SAS_ExAC"] = '0.305'
        maxFreqString = selectMaxAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(maxFreqString, '-')

        for attr, value in self.newRowAlleleFrequencies.iteritems():
            self.newRowAlleleFrequencies[attr] = '-'

        maxFreqString = selectMaxAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(maxFreqString, '-')

    def test_select_allele_frequency(self):
        AF = selectAlleleFrequency(self.newRowAlleleFrequencies)
        self.assertEquals(AF, '0.20345 (ESP)')


if __name__ == '__main__':
    pass
