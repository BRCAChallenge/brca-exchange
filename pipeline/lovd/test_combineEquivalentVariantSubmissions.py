import pytest
import unittest
from combineEquivalentVariantSubmissions import mergeRows


class TestStringMethods(unittest.TestCase):

    def setUp(self):
        self.oldRow = {
            'submission_id': 'NM_007294.3:c.(5467+1_5468-1)_(*1_?)delJohan den Dunnen (Rotterdam,NL)',
            'individuals': '1',
            'genetic_origin': 'Germline',
            'DBID': 'BRCA1_001694',
            'RNA': 'r.?',
            'edited_date': '2018-04-02 11:41:21',
            'gDNA': 'g.(?_41197694)_(41197820_41199659)del',
            'frequency': '',
            'chromosome': 'chr17',
            'cDNA': 'NM_007294.3:c.(5467+1_5468-1)_(*1_?)del',
            'created_date': '2018-04-02 11:41:21',
            'Protein': 'p.?',
            'submitters': 'Johan den Dunnen (Rotterdam,NL)',
            'geneid': 'BRCA1',
            'variant_effect': '+/.',
            'remarks': ''
        }
        self.newRow = {
            'submission_id': 'NM_007294.3:c.(5467+1_5468-1)_(*1_?)delJohan den Dunnen (Rotterdam,NL)',
            'individuals': '1',
            'genetic_origin': 'Germline',
            'DBID': 'BRCA1_001694',
            'RNA': 'r.?',
            'edited_date': '2018-04-02 11:41:21',
            'gDNA': 'g.(?_41197694)_(41197820_41199659)del',
            'frequency': '',
            'chromosome': 'chr17',
            'cDNA': 'NM_007294.3:c.(5467+1_5468-1)_(*1_?)del',
            'created_date': '2018-04-02 11:41:21',
            'Protein': 'p.?',
            'submitters': 'Johan den Dunnen (Rotterdam,NL)',
            'geneid': 'BRCA1',
            'variant_effect': '+/.',
            'remarks': ''
        }

    def test_dont_merge_submissions_with_different_submission_ids(self):
        self.oldRow['submission_id'] = "this"
        self.newRow['submission_id'] = "that"
        self.assertRaises(AssertionError, mergeRows, self.oldRow, self.newRow)

    def test_adds_individuals_field(self):
        self.oldRow['submission_id'] = "same"
        self.oldRow['individuals'] = "24"

        self.newRow['submission_id'] = "same"
        self.newRow['individuals'] = "24"

        merged = mergeRows(self.oldRow, self.newRow)
        self.assertEqual(merged['individuals'], "48")


    def test_merges_different_values_for_same_field(self):
        self.oldRow['submission_id'] = "same"
        self.oldRow['edited_date'] = "2018-04-02 11:41:21"

        self.newRow['submission_id'] = "same"
        self.newRow['edited_date'] = "2016-06-08 12:11:18"

        merged = mergeRows(self.oldRow, self.newRow)
        self.assertIn(self.oldRow['edited_date'], merged['edited_date'])
        self.assertIn(self.newRow['edited_date'], merged['edited_date'])
        self.assertIsInstance(merged['edited_date'], list)
        self.assertEqual(len(merged['edited_date']), 2)

    def test_ignores_same_values_for_same_field(self):
        self.oldRow['submission_id'] = "same"
        self.oldRow['edited_date'] = "2018-04-02 11:41:21"

        self.newRow['submission_id'] = "same"
        self.newRow['edited_date'] = "2018-04-02 11:41:21"

        merged = mergeRows(self.oldRow, self.newRow)
        self.assertIn(self.oldRow['edited_date'], merged['edited_date'])
        self.assertIn(self.newRow['edited_date'], merged['edited_date'])
        self.assertEqual(merged['edited_date'], "2018-04-02 11:41:21")

    def test_handles_reordered_lists_correctly(self):
        self.oldRow['submission_id'] = "same"
        self.oldRow['edited_date'] = ["1", "2", "3"]

        self.newRow['submission_id'] = "same"
        self.newRow['edited_date'] = ["3", "5", "2"]

        merged = mergeRows(self.oldRow, self.newRow)
        self.assertEqual(len(merged['edited_date']), 4)
        self.assertIn("1", merged['edited_date'])
        self.assertIn("2", merged['edited_date'])
        self.assertIn("3", merged['edited_date'])
        self.assertIn("5", merged['edited_date'])

    def test_only_updates_different_fields_and_individuals(self):
        oldRow = {
            'submission_id': 'same',
            'edited_date': 'this',
            'individuals': '1',
            'variant_effect': '+/.'
        }

        newRow = {
            'submission_id': 'same',
            'edited_date': 'that',
            'individuals': '1',
            'variant_effect': '+/.'
        }

        merged = mergeRows(oldRow, newRow)

        # ignores same values for variant_effect
        self.assertIsInstance(merged['variant_effect'], basestring)
        self.assertEqual(merged['variant_effect'], '+/.')

        # adds individuals
        self.assertEqual(merged['individuals'], '2')

        # merges edited_date
        self.assertIsInstance(merged['edited_date'], list)
        self.assertEqual(len(merged['edited_date']), 2)
        self.assertIn(oldRow['edited_date'], merged['edited_date'])
        self.assertIn(newRow['edited_date'], merged['edited_date'])

if __name__ == '__main__':
    pass
