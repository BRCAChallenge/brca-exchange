import pytest
import unittest
from releaseDiff import equivalentPathogenicityAllValues, checkPathogenicityAllDiffBySource
import pdb


class TestStringMethods(unittest.TestCase):

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

    def test_pathogenicity_all_diff_by_source_same_values_different_order(self):
        prev = "Uncertain_significance,Likely_benign (ClinVar)"
        new = "Likely_benign,Uncertain_significance (ClinVar)"
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertIs(classificationAdded, '')
        self.assertIs(classificationRemoved, '')

    def test_pathogenicity_all_diff_by_source_change(self):
        prev = ["Uncertain_significance,Likely_benign (ClinVar)", "Pending (BIC)"]
        new = ["Likely_benign,Pathogenic (ClinVar)", "Pending (BIC)"]

        # Test ClinVar change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("ClinVar", prev, new)
        self.assertEqual(classificationAdded, 'Pathogenic (ClinVar)')
        self.assertEqual(classificationRemoved, 'Uncertain_significance (ClinVar)')

        # Test BIC no change
        (classificationAdded, classificationRemoved) = checkPathogenicityAllDiffBySource("BIC", prev, new)
        self.assertIs(classificationAdded, '')
        self.assertIs(classificationRemoved, '')

if __name__ == '__main__':
    pass
