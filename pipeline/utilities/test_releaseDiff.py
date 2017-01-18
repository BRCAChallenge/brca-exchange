import pytest
import unittest
from releaseDiff import equivalentPathogenicityAllValues, checkPathogenicityAllDiffBySource, determineDiffForPathogenicityAll
import pdb


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

if __name__ == '__main__':
    pass
