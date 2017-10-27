import pytest
import unittest
import tempfile
import csv
import releaseDiff
import copy
from os import path


class TestStringMethods(unittest.TestCase):

    def setUp(self):
        self.fieldnames = [
                      'Pathogenicity_all',
                      'Pathogenicity_expert',
                      'Source',
                      'HGVS_Protein',
                      'Pathogenicity_expert',
                      'Genomic_Coordinate_hg38',
                      'pyhgvs_Genomic_Coordinate_38',
                      'pyhgvs_Genomic_Coordinate_37',
                      'pyhgvs_Genomic_Coordinate_36',
                      'pyhgvs_Protein',
                      'Submitter_ClinVar',
                      'Source_URL',
                      'Clinical_Significance_ClinVar',
                      'Allele_frequency_ExAC',
                      'EAS_Allele_frequency_1000_Genomes',
                      'Polyphen_Prediction',
                      'Polyphen_Score',
                      'Allele_count_AFR',
                      'Allele_Frequency',
                      'Allele_frequency_FIN_ExAC',
                      'Max_Allele_Frequency'
                     ]
        self.oldRow = {
                  'Pathogenicity_all': '',
                  'Pathogenicity_expert': '',
                  'Source': 'LOVD,1000_Genomes',
                  'HGVS_Protein': '-',
                  'Pathogenicity_expert': 'Not Yet Reviewed',
                  'Genomic_Coordinate_hg38': 'chr17:g.43049067:C>T',
                  'pyhgvs_Genomic_Coordinate_38': 'chr17:g.43049067:C>T',
                  'pyhgvs_Genomic_Coordinate_37': 'chr17:g.43049067:C>T',
                  'pyhgvs_Genomic_Coordinate_36': 'chr17:g.43049067:C>T',
                  'pyhgvs_Protein': 'NP_009225.1:p.?',
                  'Allele_frequency_ExAC': '9.841E-06',
                  'EAS_Allele_frequency_1000_Genomes': '0',
                  'Polyphen_Prediction': 'benign',
                  'Polyphen_Score': '0.992',
                  'Allele_count_AFR': '-',
                  'Allele_Frequency': '-',
                  'Allele_frequency_FIN_ExAC': '-',
                  'Max_Allele_Frequency': '0.000132 (EAS from ExAC minus TCGA)'
                 }

        self.newRow = copy.deepcopy(self.oldRow)

        self.test_dir = tempfile.mkdtemp()

        self.added = csv.DictWriter(open(path.join(self.test_dir, 'added.tsv'), 'w'), delimiter="\t", fieldnames=self.fieldnames)
        self.added.writeheader()

        self.removed = csv.DictWriter(open(path.join(self.test_dir, 'removed.tsv'), 'w'), delimiter="\t", fieldnames=self.fieldnames)
        self.removed.writeheader()

        self.added_data = open(path.join(self.test_dir, 'added_data.tsv'), 'w')

        self.diff = open(path.join(self.test_dir, 'diff.txt'), 'w')

        self.total_variants_with_additions = 0

        self.total_variants_with_changes = 0

        self.diff_json = {}

    ###################################
    # Tests for determining change type
    ###################################

    def test_reordered_pathogenicity_all_data(self):
        prev = "Pathogenic(ENIGMA); not_provided,Pathogenic (ClinVar); Class 5 (BIC)"
        new = "Pathogenic(ENIGMA); Pathogenic,not_provided (ClinVar); Class 5 (BIC)"
        prevTwo = "Uncertain_significance,Likely_benign (ClinVar); Pending (BIC)"
        newTwo = "Likely_benign,Uncertain_significance (ClinVar); Pending (BIC)"
        added, removed = releaseDiff.determineDiffForPathogenicityAll(prev, new)
        addedTwo, removedTwo = releaseDiff.determineDiffForPathogenicityAll(prevTwo, newTwo)
        self.assertIsNone(added)
        self.assertIsNone(removed)
        self.assertIsNone(addedTwo)
        self.assertIsNone(addedTwo)

    def test_swapped_pathogenicity_all_data(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar); Pending (BIC)"
        new = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        added, removed = releaseDiff.determineDiffForPathogenicityAll(prev, new)
        self.assertIn('Uncertain_significance,Likely_benign (BIC)', added)
        self.assertIn('Pending (ClinVar)', added)
        self.assertIn('Pending (BIC)', removed)
        self.assertIn('Uncertain_significance,Likely_benign (ClinVar)', removed)

    def test_different_pathogenicity_all_data(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar); Pending (BIC)"
        new = "Likely_benign (ClinVar); Pending (BIC)"
        added, removed = releaseDiff.determineDiffForPathogenicityAll(prev, new)
        self.assertIsNone(added)
        self.assertIn('Uncertain_significance (ClinVar)', removed)

    def test_same_pathogenicity_all_data_single_source(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar)"
        new = "Likely_benign,Uncertain_significance (ClinVar)"
        added, removed = releaseDiff.determineDiffForPathogenicityAll(prev, new)
        self.assertIsNone(added)
        self.assertIsNone(removed)

    ###################################
    # Tests for determining diff json
    ###################################

    # Pathogenicity_all Tests:

    def test_pathogenicity_all_diff_by_source_same_values_different_order(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar)"
        new = "Likely_benign,Uncertain_significance (ClinVar)"
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, '')
        self.assertEqual(classificationRemoved, '')

    def test_pathogenicity_all_diff_by_source_change(self):
        prev = ["Uncertain_significance,Likely_benign (ClinVar)", "Pending (BIC)"]
        new = ["Likely_benign,Pathogenic (ClinVar)", "Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, 'Pathogenic (ClinVar)')
        self.assertEqual(classificationRemoved, 'Uncertain_significance (ClinVar)')

        # Test BIC no change
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertEqual(classificationAdded, '')
        self.assertEqual(classificationRemoved, '')

    def test_pathogenicity_all_diff_by_source_swap_sources_changed_classifications(self):
        prev = ["Uncertain_significance,Likely_benign (BIC)", "Pending (ClinVar)"]
        new = ["Uncertain_significance,Likely_benign,Pending (ClinVar)", "Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, 'Uncertain_significance,Likely_benign (ClinVar)')
        self.assertEqual(classificationRemoved, '')

        # Test BIC
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertEqual(classificationAdded, 'Pending (BIC)')
        self.assertEqual(classificationRemoved, 'Uncertain_significance,Likely_benign (BIC)')

    def test_pathogenicity_all_diff_by_source_drop_source_add_source(self):
        prev = ["Pending (ClinVar)"]
        new = ["Uncertain_significance,Likely_benign,Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, '')
        self.assertEqual(classificationRemoved, 'Pending (ClinVar)')

        # Test BIC change
        (classificationAdded, classificationRemoved) = releaseDiff.checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertEqual(classificationAdded, 'Uncertain_significance,Likely_benign,Pending (BIC)')
        self.assertEqual(classificationRemoved, '')

    def test_pathogenicity_all_diff_change(self):
        prev = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        new = "Uncertain_significance,Likely_benign,Pending (ClinVar); Pending (BIC)"
        (added, removed) = releaseDiff.determineDiffForPathogenicityAll(prev, new)
        self.assertIn('Pending (BIC)', added)
        self.assertIn('Uncertain_significance,Likely_benign (ClinVar)', added)
        self.assertIn('Uncertain_significance,Likely_benign (BIC)', removed)
        self.assertNotIn('Pending (ClinVar)', removed)

    def test_pathogenicity_all_diff_no_change(self):
        prev = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        new = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        (added, removed) = releaseDiff.determineDiffForPathogenicityAll(prev, new)
        self.assertIsNone(added)
        self.assertIsNone(removed)

    def test_diff_json_handles_no_changes_correctly(self):
        field = "Pathogenicity_all"
        prev = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        new = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Pathogenicity_all')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIsNone(diff['added'])
        self.assertIsNone(diff['removed'])

    # General diff tests:

    def test_diff_json_handles_list_changes_correctly(self):
        field = "Source"
        prev = "BIC,ClinVar,ENIGMA"
        new = "ClinVar"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Source')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIsNone(diff['added'])
        self.assertEqual(diff['removed'], ['BIC', 'ENIGMA'])

    def test_diff_json_handles_list_changes_correctly_changed_order(self):
        field = "Source_URL"
        prev = "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000075545, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000187590, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000144137"
        new = "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000187590, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000075545, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000144137"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Source_URL')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIsNone(diff['added'])
        self.assertIsNone(diff['removed'])

    def test_diff_json_handles_list_data_added_from_nothing_correctly(self):
        field = "Submitter_ClinVar"
        prev = "-"
        new = ("c/o_University_of_Cambridge,The_Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA),"
               "Ambry_Genetics,Breast_Cancer_Information_Core_(BIC)_(BRCA2),Invitae")
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Submitter_ClinVar')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIn('c/o_University_of_Cambridge', diff['added'])
        self.assertIn('The_Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA)', diff['added'])
        self.assertIn('Ambry_Genetics', diff['added'])
        self.assertIn('Breast_Cancer_Information_Core_(BIC)_(BRCA2)', diff['added'])
        self.assertIn('Invitae', diff['added'])
        self.assertEqual(diff['removed'], ['-'])

    def test_diff_json_handles_individual_changes_correctly(self):
        field = "HGVS_Protein"
        prev = "NM_000059:p.His1085Arg"
        new = "p.(His1085Arg)"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertEqual(diff['added'], 'p.(His1085Arg)')
        self.assertEqual(diff['removed'], 'NM_000059:p.His1085Arg')

    def test_diff_json_handles_individual_changes_data_removed_correctly(self):
        field = "HGVS_Protein"
        prev = "NM_000059:p.His1085Arg"
        new = "-"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertEqual(diff['added'], '-')
        self.assertEqual(diff['removed'], 'NM_000059:p.His1085Arg')

    def test_diff_json_handles_individual_changes_data_added_from_nothing_correctly(self):
        field = "HGVS_Protein"
        prev = "-"
        new = "NM_000059:p.His1085Arg"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertEqual(diff['added'], 'NM_000059:p.His1085Arg')
        self.assertEqual(diff['removed'], '-')

    def test_diff_json_handles_hgvs_protein_non_changes_correctly(self):
        field = "HGVS_Protein"
        prev = "NM_000059:p.His1085Arg"
        new = "NM_000059:p.His1085Arg"
        diff = releaseDiff.determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertIsNone(diff['added'])
        self.assertIsNone(diff['removed'])

    # Test full releaseDiff output by using rows of data from release

    def test_compare_row_equal_rows(self):
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        self.assertIsNone(change_type)

    def test_compare_row_added_data(self):
        releaseDiff.added_data = self.added_data
        self.newRow['Source'] += ",ENIGMA"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        self.assertEqual(change_type, "added_information")

    def test_ignores_unused_column_changes_in_compare_row(self):
        releaseDiff.added_data = self.added_data
        self.oldRow['pyhgvs_Protein'] = "NP_009225.1:p.?"
        self.newRow['HGVS_Protein'] = "NM_000059:p.His1085Arg"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        self.assertIsNone(change_type)

    def test_catches_new_column_changes_in_compare_row(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow['pyhgvs_Protein'] = "NP_009225.1:p.?"
        self.newRow['pyhgvs_Protein'] = "NM_000059:p.His1085Arg"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json[variant]
        self.assertEqual(change_type, "changed_information")
        self.assertIs(diff[0]['removed'], 'NP_009225.1:p.?')
        self.assertIs(diff[0]['added'], 'NM_000059:p.His1085Arg')

    def test_builds_diff_with_adjusted_column_names_to_match_db(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow['pyhgvs_Protein'] = "NP_009225.1:p.?"
        self.newRow['pyhgvs_Protein'] = "NM_000059:p.His1085Arg"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json[variant]
        self.assertIs(diff[0]['field'], 'HGVS_Protein')

    def test_catches_pathogenicity_changes(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow['Pathogenicity_all'] = "Pathogenic,not_provided (ClinVar); Class 5 (BIC)"
        self.oldRow['Pathogenicity_expert'] = "Not Yet Classified"
        self.newRow['Pathogenicity_all'] = "Pathogenic(ENIGMA); Pathogenic,not_provided (ClinVar); Class 5 (BIC)"
        self.newRow['Pathogenicity_expert'] = "Pathogenic"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json[variant]
        self.assertEqual(len(diff), 2)
        self.assertIs(change_type, "changed_classification")

    def test_change_type_ignores_reorders(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow['Pathogenicity_all'] = "Pathogenic(ENIGMA); not_provided,Pathogenic (ClinVar); Class 5 (BIC)"
        self.oldRow['Pathogenicity_expert'] = "Pathogenic"
        self.newRow['Pathogenicity_all'] = "Pathogenic(ENIGMA); Pathogenic,not_provided (ClinVar); Class 5 (BIC)"
        self.newRow['Pathogenicity_expert'] = "Pathogenic"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(len(diff), 0)
        self.assertIsNone(change_type)

    def test_add_gs_to_genomic_coordinate(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        self.oldRow['pyhgvs_Genomic_Coordinate_38'] = "chr17:43049067:C>T"
        self.newRow['pyhgvs_Genomic_Coordinate_38'] = "chr17:g.43049067:C>T"
        self.oldRow['pyhgvs_Genomic_Coordinate_37'] = "chr17:43049067:C>T"
        self.newRow['pyhgvs_Genomic_Coordinate_37'] = "chr17:g.43049067:C>T"
        self.oldRow['pyhgvs_Genomic_Coordinate_36'] = "chr17:43049067:C>T"
        self.newRow['pyhgvs_Genomic_Coordinate_36'] = "chr17:g.43049067:C>T"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        self.oldRow = releaseDiff.addGsIfNecessary(self.oldRow)
        self.oldRow = releaseDiff.addGsIfNecessary(self.newRow)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        self.assertIsNone(change_type)

    def test_properly_classifies_variants_with_removed_columns(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow["Functional_analysis_result_LOVD"] = "test"
        self.newRow['Source'] += ",ENIGMA"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json[variant]
        self.assertEqual(len(diff), 2)
        self.assertIs(change_type, "changed_information")

    def test_properly_classifies_variants_with_removed_columns_of_empty_data(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow["Functional_analysis_result_LOVD"] = "-"
        self.newRow['Source'] += ",ENIGMA"
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json[variant]
        self.assertEqual(len(diff), 1)
        self.assertIs(change_type, "added_information")

    def test_ignores_cosmetic_changes_in_diff(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)

        self.oldRow["Submitter_ClinVar"] = "The_Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA),c/o_University_of_Cambridge"
        self.newRow["Submitter_ClinVar"] = "Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA),_c/o_University_of_Cambridge"

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

        self.oldRow["Submitter_ClinVar"] = "Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA),_c/o_University_of_Cambridge"
        self.newRow["Submitter_ClinVar"] = "The_Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA),c/o_University_of_Cambridge"

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

        self.newRow['EAS_Allele_frequency_1000_Genomes'] = '0.0'
        self.newRow['Allele_frequency_ExAC'] = '9.841e-06'

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

        self.oldRow['EAS_Allele_frequency_1000_Genomes'] = '0.0'
        self.oldRow['Allele_frequency_ExAC'] = '9.841e-06'
        self.newRow['EAS_Allele_frequency_1000_Genomes'] = '0'
        self.newRow['Allele_frequency_ExAC'] = '9.841E-06'

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_catches_changed_numeric_values_after_normalization(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)

        self.oldRow['Allele_frequency_ExAC'] = '9.841e-06'
        self.newRow['Allele_frequency_ExAC'] = '9.841e-07'

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(len(diff), 1)
        self.assertIs(change_type, "changed_information")

    def test_catches_reordered_source_urls(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        self.oldRow["Source_URL"] = "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000075538, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000144133, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000109288"
        self.newRow["Source_URL"] = "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000144133, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000075538, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000109288"
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_handle_repeat_data_correctly(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        self.oldRow["Pathogenicity_all"] = "Pathogenic,Pathogenic,not_provided"
        self.newRow["Pathogenicity_all"] = "Pathogenic,not_provided"
        self.oldRow["Clinical_Significance_ClinVar"] = "Pathogenic,Pathogenic,not_provided"
        self.newRow["Clinical_Significance_ClinVar"] = "Pathogenic,not_provided"
        self.oldRow["Submitter_ClinVar"] = "PreventionGenetics"
        self.newRow["Submitter_ClinVar"] = "PreventionGenetics,PreventionGenetics"
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_handles_new_data_in_new_columns_in_v2_data_correctly(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.updated_fieldnames = self.fieldnames + ['Genetic_origin_LOVD', 'RNA_LOVD', 'Submitters_LOVD']
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.updated_fieldnames)
        self.newRow["Genetic_origin_LOVD"] = "lorem ipsum"
        self.newRow["RNA_LOVD"] = "lorem ipsum"
        self.newRow["Submitters_LOVD"] = "lorem ipsum"
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(change_type, "added_information")
        self.assertTrue(len(diff[variant]) == 3)

    def test_handles_null_data_in_new_columns_in_v2_data_correctly(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.updated_fieldnames = self.fieldnames + ['Genetic_origin_LOVD', 'RNA_LOVD', 'Submitters_LOVD']
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.updated_fieldnames)
        self.newRow["Genetic_origin_LOVD"] = ""
        self.newRow["RNA_LOVD"] = None
        self.newRow["Submitters_LOVD"] = "-"
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_ignores_polyphen_fields(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)

        self.newRow["Polyphen_Score"] = "0.283"
        self.newRow["Polyphen_Prediction"] = "probably_damaging"

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_ignores_max_allele_frequency_field(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)

        self.newRow["Max_Allele_Frequency"] = "-"

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_ignores_change_from_none_to_empty_string(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)

        self.oldRow["Submitter_ClinVar"] = "Quest_Diagnostics_Nichols_Institute_San_Juan_Capistrano"
        self.newRow["Submitter_ClinVar"] = ",Quest_Diagnostics_Nichols_Institute_San_Juan_Capistrano"

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_handles_fields_with_commas_in_parentheses_correctly(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.updated_fieldnames = self.fieldnames + ['Genetic_origin_LOVD', 'RNA_LOVD', 'Submitters_LOVD']
        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.updated_fieldnames)

        self.newRow["Submitters_LOVD"] = "Ans M.W. van den Ouweland (Rotterdam,NL), Genevieve Michils (Leuven,BE), Rien Blok (Maastricht NL)"

        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(len(diff), 1)
        v_diff = diff['chr17:g.43049067:C>T'][0]
        self.assertEqual(v_diff['field'], 'Submitters_LOVD')
        self.assertEqual(len(v_diff['added']), 3)
        self.assertIn("Ans M.W. van den Ouweland (Rotterdam,NL)", v_diff['added'])
        self.assertIn("Genevieve Michils (Leuven,BE)", v_diff['added'])
        self.assertIn("Rien Blok (Maastricht NL)", v_diff['added'])

    def test_displays_integers_as_integers_not_floats(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'
        self.oldRow['Allele_count_AFR'] = '-'
        self.newRow['Allele_count_AFR'] = '567'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        v_diff = diff['chr17:g.43049067:C>T'][0]
        self.assertEqual(len(diff), 1)
        self.assertEqual(v_diff['field'], 'Allele_count_AFR')
        self.assertEqual(v_diff['added'], '567')
        self.assertEqual(v_diff['removed'], '-')

    def test_properly_handles_field_name_changes_same_data(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.updated_fieldnames = self.fieldnames + ['Minor_allele_frequency_percent_ESP']
        self.fieldnames = self.fieldnames + ['Minor_allele_frequency_ESP']

        self.oldRow['Minor_allele_frequency_ESP'] = '2.5'
        self.newRow['Minor_allele_frequency_percent_ESP'] = '2.5'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.updated_fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_properly_handles_field_name_changes_updated_data(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.updated_fieldnames = self.fieldnames + ['Minor_allele_frequency_percent_ESP']
        self.fieldnames = self.fieldnames + ['Minor_allele_frequency_ESP']

        self.oldRow['Minor_allele_frequency_ESP'] = '1.5'
        self.newRow['Minor_allele_frequency_percent_ESP'] = '2.5'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.updated_fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        v_diff = diff['chr17:g.43049067:C>T'][0]
        self.assertEqual(len(diff), 1)
        self.assertEqual(v_diff['field'], 'Minor_allele_frequency_percent_ESP')
        self.assertEqual(v_diff['added'], '2.5')
        self.assertEqual(v_diff['removed'], '1.5')
        self.assertEqual(change_type, "changed_information")

    def test_properly_handles_field_name_changes_added_data(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.updated_fieldnames = self.fieldnames + ['Minor_allele_frequency_percent_ESP']
        self.fieldnames = self.fieldnames + ['Minor_allele_frequency_ESP']

        self.oldRow['Minor_allele_frequency_ESP'] = '-'
        self.newRow['Minor_allele_frequency_percent_ESP'] = '2.5'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.updated_fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        v_diff = diff['chr17:g.43049067:C>T'][0]
        self.assertEqual(len(diff), 1)
        self.assertEqual(v_diff['field'], 'Minor_allele_frequency_percent_ESP')
        self.assertEqual(v_diff['added'], '2.5')
        self.assertEqual(v_diff['removed'], '-')
        self.assertEqual(change_type, "added_information")

    def test_ignored_exac_minus_tcga_instead_of_exac(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.oldRow['Allele_Frequency'] = '9.42e-06 (ExAC)'
        self.newRow['Allele_Frequency'] = '9.42e-06 (ExAC minus TCGA)'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_catches_value_change_for_exac_allele_frequency(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.oldRow['Allele_Frequency'] = '9.42e-06 (ExAC)'
        self.newRow['Allele_Frequency'] = '9.99e-06 (ExAC minus TCGA)'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        v_diff = diff['chr17:g.43049067:C>T'][0]
        self.assertEqual(len(diff), 1)
        self.assertEqual(v_diff['field'], 'Allele_Frequency')
        self.assertEqual(v_diff['added'], '9.99e-06 (ExAC minus TCGA)')
        self.assertEqual(v_diff['removed'], '9.42e-06 (ExAC minus TCGA)')
        self.assertEqual(change_type, "changed_information")

    def test_properly_ignores_exac_minus_tcga_rounding_changes_for_generic_allele_frequency(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.oldRow['Allele_Frequency'] = '9.419e-06 (ExAC minus TCGA)'
        self.newRow['Allele_Frequency'] = '9.42e-06 (ExAC minus TCGA)'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_ignores_rounding_for_exac_af_fields(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.oldRow['Allele_frequency_FIN_ExAC'] = '9.4181237657'
        self.newRow['Allele_frequency_FIN_ExAC'] = '9.42'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        self.assertEqual(diff, {})
        self.assertIsNone(change_type)

    def test_different_values_for_exac_af_fields(self):
        releaseDiff.added_data = self.added_data
        releaseDiff.diff = self.diff
        releaseDiff.diff_json = self.diff_json
        variant = 'chr17:g.43049067:C>T'

        self.oldRow['Allele_frequency_FIN_ExAC'] = '9.4181237657'
        self.newRow['Allele_frequency_FIN_ExAC'] = '9.41'

        v1v2 = releaseDiff.v1ToV2(self.fieldnames, self.fieldnames)
        change_type = v1v2.compareRow(self.oldRow, self.newRow)
        diff = releaseDiff.diff_json
        v_diff = diff['chr17:g.43049067:C>T'][0]
        self.assertEqual(len(diff), 1)
        self.assertEqual(v_diff['field'], 'Allele_frequency_FIN_ExAC')
        self.assertEqual(v_diff['added'], '9.41')
        self.assertEqual(v_diff['removed'], '9.42')
        self.assertEqual(change_type, "changed_information")


if __name__ == '__main__':
    pass
