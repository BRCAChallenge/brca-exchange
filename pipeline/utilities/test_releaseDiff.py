import pytest
import unittest
from releaseDiff import (equivalentPathogenicityAllValues, checkPathogenicityAllDiffBySource,
                         determineDiffForPathogenicityAll, determineDiffForJSON, transformer,
                         v1ToV2)


class TestStringMethods(unittest.TestCase):

    ###################################
    # Tests for determining change type
    ###################################

    def test_reordered_pathogenicity_all_data(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar); Pending (BIC)"
        new = "Likely_benign,Uncertain_significance (ClinVar); Pending (BIC)"
        self.assertTrue(equivalentPathogenicityAllValues(prev, new))

    def test_swapped_pathogenicity_all_data(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar); Pending (BIC)"
        new = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        self.assertFalse(equivalentPathogenicityAllValues(prev, new))

    def test_different_pathogenicity_all_data(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar); Pending (BIC)"
        new = "Likely_benign (ClinVar); Pending (BIC)"
        self.assertFalse(equivalentPathogenicityAllValues(prev, new))

    def test_same_pathogenicity_all_data_single_source(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar)"
        new = "Likely_benign,Uncertain_significance (ClinVar)"
        self.assertTrue(equivalentPathogenicityAllValues(prev, new))

    ###################################
    # Tests for determining diff json
    ###################################

    # Pathogenicity_all Tests:

    def test_pathogenicity_all_diff_by_source_same_values_different_order(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar)"
        new = "Likely_benign,Uncertain_significance (ClinVar)"
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, '')
        self.assertEqual(classificationRemoved, '')

    def test_pathogenicity_all_diff_by_source_change(self):
        prev = ["Uncertain_significance,Likely_benign (ClinVar)", "Pending (BIC)"]
        new = ["Likely_benign,Pathogenic (ClinVar)", "Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, 'Pathogenic (ClinVar)')
        self.assertEqual(classificationRemoved, 'Uncertain_significance (ClinVar)')

        # Test BIC no change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertEqual(classificationAdded, '')
        self.assertEqual(classificationRemoved, '')

    def test_pathogenicity_all_diff_by_source_swap_sources_changed_classifications(self):
        prev = ["Uncertain_significance,Likely_benign (BIC)", "Pending (ClinVar)"]
        new = ["Uncertain_significance,Likely_benign,Pending (ClinVar)", "Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, 'Uncertain_significance,Likely_benign (ClinVar)')
        self.assertEqual(classificationRemoved, '')

        # Test BIC
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertEqual(classificationAdded, 'Pending (BIC)')
        self.assertEqual(classificationRemoved, 'Uncertain_significance,Likely_benign (BIC)')

    def test_pathogenicity_all_diff_by_source_drop_source_add_source(self):
        prev = ["Pending (ClinVar)"]
        new = ["Uncertain_significance,Likely_benign,Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, '')
        self.assertEqual(classificationRemoved, 'Pending (ClinVar)')

        # Test BIC change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertEqual(classificationAdded, 'Uncertain_significance,Likely_benign,Pending (BIC)')
        self.assertEqual(classificationRemoved, '')

    def test_pathogenicity_all_diff_change(self):
        prev = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        new = "Uncertain_significance,Likely_benign,Pending (ClinVar); Pending (BIC)"
        (added, removed) = determineDiffForPathogenicityAll(prev, new)
        self.assertIn('Pending (BIC)', added)
        self.assertIn('Uncertain_significance,Likely_benign (ClinVar)', added)
        self.assertIn('Uncertain_significance,Likely_benign (BIC)', removed)
        self.assertNotIn('Pending (ClinVar)', removed)

    def test_pathogenicity_all_diff_no_change(self):
        prev = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        new = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        (added, removed) = determineDiffForPathogenicityAll(prev, new)
        self.assertIsNone(added)
        self.assertIsNone(removed)

    def test_diff_json_handles_no_changes_correctly(self):
        field = "Pathogenicity_all"
        prev = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        new = "Uncertain_significance,Likely_benign (BIC); Pending (ClinVar)"
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Pathogenicity_all')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIsNone(diff['added'])
        self.assertIsNone(diff['removed'])

    # General diff tests:

    def test_diff_json_handles_list_changes_correctly(self):
        field = "Source"
        prev = "BIC,ClinVar,ENIGMA"
        new = "ClinVar"
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Source')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIsNone(diff['added'])
        self.assertEqual(diff['removed'], ['BIC', 'ENIGMA'])

    def test_diff_json_handles_list_changes_correctly_changed_order(self):
        field = "Source_URL"
        prev = "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000075545, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000187590, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000144137"
        new = "http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000187590, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000075545, http://www.ncbi.nlm.nih.gov/clinvar/?term=SCV000144137"
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'Source_URL')
        self.assertEqual(diff['field_type'], 'list')
        self.assertIsNone(diff['added'])
        self.assertIsNone(diff['removed'])

    def test_diff_json_handles_list_data_added_from_nothing_correctly(self):
        field = "Submitter_ClinVar"
        prev = "-"
        new = ("c/o_University_of_Cambridge,The_Consortium_of_Investigators_of_Modifiers_of_BRCA1/2_(CIMBA),"
               "Ambry_Genetics,Breast_Cancer_Information_Core_(BIC)_(BRCA2),Invitae")
        diff = determineDiffForJSON(field, prev, new)
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
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertEqual(diff['added'], 'p.(His1085Arg)')
        self.assertEqual(diff['removed'], 'NM_000059:p.His1085Arg')

    def test_diff_json_handles_individual_changes_data_removed_correctly(self):
        field = "HGVS_Protein"
        prev = "NM_000059:p.His1085Arg"
        new = "-"
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertEqual(diff['added'], '-')
        self.assertEqual(diff['removed'], 'NM_000059:p.His1085Arg')

    def test_diff_json_handles_individual_changes_data_added_from_nothing_correctly(self):
        field = "HGVS_Protein"
        prev = "-"
        new = "NM_000059:p.His1085Arg"
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertEqual(diff['added'], 'NM_000059:p.His1085Arg')
        self.assertEqual(diff['removed'], '-')

    def test_diff_json_handles_hgvs_protein_non_changes_correctly(self):
        field = "HGVS_Protein"
        prev = "NM_000059:p.His1085Arg"
        new = "NM_000059:p.His1085Arg"
        diff = determineDiffForJSON(field, prev, new)
        self.assertEqual(diff['field'], 'HGVS_Protein')
        self.assertEqual(diff['field_type'], 'individual')
        self.assertIsNone(diff['added'])
        self.assertIsNone(diff['removed'])

    def test_compare_row_equal_rows(self):
        fieldnames = [
                      'Pathogenicity_all',
                      'pyhgvs_cDNA',
                      'Source',
                      'HGVS_Protein',
                      'Polyphen_Score',
                      'pyhgvs_Hg37_Start',
                      'Hg38_End',
                      'HGVS_cDNA',
                      'Hg37_End',
                      'Hg38_Start',
                      'Pathogenicity_expert',
                      'pyhgvs_Hg36_End',
                      'Hg36_End',
                      'pyhgvs_Hg37_End',
                      'Genomic_Coordinate_hg37',
                      'Genomic_Coordinate_hg36',
                      'pyhgvs_Genomic_Coordinate_37',
                      'pyhgvs_Genomic_Coordinate_36',
                      'Genomic_Coordinate_hg38',
                      'pyhgvs_Hg36_Start',
                      'pyhgvs_Genomic_Coordinate_38',
                      'Protein_Change',
                      'HGVS_RNA',
                      'Hg37_Start',
                      'Hg36_Start',
                      'pyhgvs_Protein'
                     ]
        oldRow = {
                  'Pathogenicity_all': '',
                  'pyhgvs_cDNA': 'NM_007294.3:c.5406+54G>A',
                  'Source': 'LOVD,1000_Genomes',
                  'HGVS_Protein': '-',
                  'Polyphen_Score': '-',
                  'pyhgvs_Hg37_Start': '41201084',
                  'Hg38_End': '43049067',
                  'HGVS_cDNA': 'c.5406+54G>A',
                  'Hg37_End': '-',
                  'Hg38_Start': '43049067',
                  'Pathogenicity_expert': 'Not Yet Reviewed',
                  'pyhgvs_Hg36_End': '38454610',
                  'Hg36_End': '-',
                  'pyhgvs_Hg37_End': '41201084',
                  'Genomic_Coordinate_hg37': '-',
                  'Genomic_Coordinate_hg36': '-',
                  'pyhgvs_Genomic_Coordinate_37': 'chr17:g.41201084:C>T',
                  'pyhgvs_Genomic_Coordinate_36': 'chr17:g.38454610:C>T',
                  'Genomic_Coordinate_hg38': 'chr17:g.43049067:C>T',
                  'pyhgvs_Hg36_Start': '38454610',
                  'pyhgvs_Genomic_Coordinate_38': 'chr17:g.43049067:C>T',
                  'Protein_Change': '-',
                  'HGVS_RNA': '-',
                  'Hg37_Start': '-',
                  'Hg36_Start': '-',
                  'pyhgvs_Protein': 'NP_009225.1:p.?'
                 }
        newRow = {
                  'Pathogenicity_all': '',
                  'pyhgvs_cDNA': 'NM_007294.3:c.5406+54G>A',
                  'Source': 'LOVD,1000_Genomes',
                  'HGVS_Protein': '-',
                  'Polyphen_Score': '-',
                  'pyhgvs_Hg37_Start': '41201084',
                  'Hg38_End': '43049067',
                  'HGVS_cDNA': 'c.5406+54G>A',
                  'Hg37_End': '-',
                  'Hg38_Start': '43049067',
                  'Pathogenicity_expert': 'Not Yet Reviewed',
                  'pyhgvs_Hg36_End': '38454610',
                  'Hg36_End': '-',
                  'pyhgvs_Hg37_End': '41201084',
                  'Genomic_Coordinate_hg37': '-',
                  'Genomic_Coordinate_hg36': '-',
                  'pyhgvs_Genomic_Coordinate_37': 'chr17:g.41201084:C>T',
                  'pyhgvs_Genomic_Coordinate_36': 'chr17:g.38454610:C>T',
                  'Genomic_Coordinate_hg38': 'chr17:g.43049067:C>T',
                  'pyhgvs_Hg36_Start': '38454610',
                  'pyhgvs_Genomic_Coordinate_38': 'chr17:g.43049067:C>T',
                  'Protein_Change': '-',
                  'HGVS_RNA': '-',
                  'Hg37_Start': '-',
                  'Hg36_Start': '-',
                  'pyhgvs_Protein': 'NP_009225.1:p.?'
                 }
        v1v2 = v1ToV2(fieldnames, fieldnames)
        change_type = v1v2.compareRow(oldRow, newRow)
        self.assertIsNone(change_type)

    def test_compare_row_added_data(self):
        fieldnames = [
                      'Pathogenicity_all',
                      'pyhgvs_cDNA',
                      'Source',
                      'HGVS_Protein',
                      'Polyphen_Score',
                      'pyhgvs_Hg37_Start',
                      'Hg38_End',
                      'HGVS_cDNA',
                      'Hg37_End',
                      'Hg38_Start',
                      'Pathogenicity_expert',
                      'pyhgvs_Hg36_End',
                      'Hg36_End',
                      'pyhgvs_Hg37_End',
                      'Genomic_Coordinate_hg37',
                      'Genomic_Coordinate_hg36',
                      'pyhgvs_Genomic_Coordinate_37',
                      'pyhgvs_Genomic_Coordinate_36',
                      'Genomic_Coordinate_hg38',
                      'pyhgvs_Hg36_Start',
                      'pyhgvs_Genomic_Coordinate_38',
                      'Protein_Change',
                      'HGVS_RNA',
                      'Hg37_Start',
                      'Hg36_Start',
                      'pyhgvs_Protein'
                     ]
        oldRow = {
                  'Pathogenicity_all': '',
                  'pyhgvs_cDNA': 'NM_007294.3:c.5406+54G>A',
                  'Source': 'LOVD,1000_Genomes',
                  'HGVS_Protein': '-',
                  'Polyphen_Score': '-',
                  'pyhgvs_Hg37_Start': '41201084',
                  'Hg38_End': '43049067',
                  'HGVS_cDNA': 'c.5406+54G>A',
                  'Hg37_End': '-',
                  'Hg38_Start': '43049067',
                  'Pathogenicity_expert': 'Not Yet Reviewed',
                  'pyhgvs_Hg36_End': '38454610',
                  'Hg36_End': '-',
                  'pyhgvs_Hg37_End': '41201084',
                  'Genomic_Coordinate_hg37': '-',
                  'Genomic_Coordinate_hg36': '-',
                  'pyhgvs_Genomic_Coordinate_37': 'chr17:g.41201084:C>T',
                  'pyhgvs_Genomic_Coordinate_36': 'chr17:g.38454610:C>T',
                  'Genomic_Coordinate_hg38': 'chr17:g.43049067:C>T',
                  'pyhgvs_Hg36_Start': '38454610',
                  'pyhgvs_Genomic_Coordinate_38': 'chr17:g.43049067:C>T',
                  'Protein_Change': '-',
                  'HGVS_RNA': '-',
                  'Hg37_Start': '-',
                  'Hg36_Start': '-',
                  'pyhgvs_Protein': 'NP_009225.1:p.?'
                 }
        newRow = {
                  'Pathogenicity_all': 'Benign (ClinVar)',
                  'pyhgvs_cDNA': 'NM_007294.3:c.5406+54G>A',
                  'Source': 'LOVD,1000_Genomes',
                  'HGVS_Protein': '-',
                  'Polyphen_Score': '-',
                  'pyhgvs_Hg37_Start': '41201084',
                  'Hg38_End': '43049067',
                  'HGVS_cDNA': 'c.5406+54G>A',
                  'Hg37_End': '-',
                  'Hg38_Start': '43049067',
                  'Pathogenicity_expert': 'Not Yet Reviewed',
                  'pyhgvs_Hg36_End': '38454610',
                  'Hg36_End': '-',
                  'pyhgvs_Hg37_End': '41201084',
                  'Genomic_Coordinate_hg37': '-',
                  'Genomic_Coordinate_hg36': '-',
                  'pyhgvs_Genomic_Coordinate_37': 'chr17:g.41201084:C>T',
                  'pyhgvs_Genomic_Coordinate_36': 'chr17:g.38454610:C>T',
                  'Genomic_Coordinate_hg38': 'chr17:g.43049067:C>T',
                  'pyhgvs_Hg36_Start': '38454610',
                  'pyhgvs_Genomic_Coordinate_38': 'chr17:g.43049067:C>T',
                  'Protein_Change': '-',
                  'HGVS_RNA': '-',
                  'Hg37_Start': '-',
                  'Hg36_Start': '-',
                  'pyhgvs_Protein': 'NP_009225.1:p.?'
                 }
        v1v2 = v1ToV2(fieldnames, fieldnames)
        change_type = v1v2.compareRow(oldRow, newRow)
        self.assertEqual(change_type, "added_information")


if __name__ == '__main__':
    pass
