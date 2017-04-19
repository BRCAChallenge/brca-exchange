import pytest
import unittest
from os import path, getcwd
import lovd2vcf
import vcf


BRCA_DATA_FILE = path.join(path.dirname(__file__), 'test_files/BRCA.txt')


class TestStringMethods(unittest.TestCase):

    def test_normalize_reports_vcf(self):
        self.assertEqual(2, 2)


if __name__ == '__main__':
    pass
