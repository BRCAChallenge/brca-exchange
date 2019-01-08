import csv
import unittest

from calcMaxEntScanMeanStd import runMaxEntScan
from calcVarPriors import getVarData
from calc_priors import extract, compute


class test_calcMaxEndScanMeanStd(unittest.TestCase):
    def setUp(self):
        pass

    def test_MSE_perl_vs_maxentpy_single(self):
        seq = 'cagGTAAGT'  # this is a donor sequence, and it should be about 10.86
        non_perl, perl = runMaxEntScan(seq, donor=True, usePerl=False), runMaxEntScan(seq, donor=True, usePerl=True)
        self.assertEqual(non_perl, 10.86)
        self.assertEqual(perl, 10.86)

    def compare_MSE(self, compare_type, long_test=False):
        """
        Uses the priors_short (or tests/priors_long) file to test runMaxEntScan's MES values against precomputed
        references. if compare_type is 'perl' uses the perl implemention, if 'python' uses the maxentpy implementation,
        and if 'cross' compares the two against each other using the sequences in the priors file, ignoring the
        precomputed MES values.

        :param long_test: if True, uses tests/priors_long.tsv as the input sequences, otherwise tests/priors_short.tsv
        :param compare_type: 'perl' to use the perl version of MES, 'python' for maxentpy, 'cross' to compare both against each other
        :return:
        """

        if compare_type not in ['perl', 'python', 'cross']:
            raise ValueError("Invalid value %s specified for compare_type, must be 'perl', 'python', or 'cross" % compare_type)

        # computes MES on column identified by first value, compares to second
        # third element is True if donor, False otherwise
        compare_cols = [
            ('refRefDonorSeq', 'refRefDonorMES', True),
            ('altRefDonorSeq', 'altRefDonorMES', True),
            ('refDeNovoDonorSeq', 'refDeNovoDonorMES', True),
            ('altDeNovoDonorSeq', 'altDeNovoDonorMES', True),
            ('closestDonorRefSeq', 'closestDonorRefMES', True),
            ('closestDonorAltSeq', 'closestDonorAltMES', True),
            ('refRefAccSeq', 'refRefAccMES', False),
            ('altRefAccSeq', 'altRefAccMES', False),
            ('refDeNovoAccSeq', 'refDeNovoAccMES', False),
            ('altDeNovoAccSeq', 'altDeNovoAccMES', False),
            ('closestAccRefSeq', 'closestAccRefMES', False),
            ('closestAccAltSeq', 'closestAccAltMES', False)
        ]

        test_file = "tests/priors_long.tsv" if long_test else "tests/priors_short.tsv"

        # look at variants in the priors_long file and check if their MES values match up to our computations
        with open(test_file, "r") as test_fp:
            test_vars = csv.DictReader(test_fp, delimiter="\t")
            for row in test_vars:
                for (seq_col, expected_MES_col, isDonor) in compare_cols:
                    seq = row[seq_col]

                    # only test variants with valid sequences
                    if seq == '-' or seq == 'N/A':
                        continue

                    if compare_type == 'cross':
                        self.assertEqual(
                            runMaxEntScan(seq, donor=isDonor, usePerl=True),
                            runMaxEntScan(seq, donor=isDonor, usePerl=False)
                        )
                    else:
                        expected_MES = row[expected_MES_col]
                        self.assertEqual(
                            runMaxEntScan(seq, donor=isDonor, usePerl=(compare_type == 'perl')),
                            float(expected_MES)
                        )

    def test_MES_maxentpy(self):
        """
        Compares the maxentpy library's MES values against MES values in tests/priors_short.tsv
        """
        self.compare_MSE(compare_type='python')

    def test_MES_perl(self):
        """
        Compares the reference perl MES values against MES values in tests/priors_short.tsv
        """
        self.compare_MSE(compare_type='perl')

    def test_MES_cross(self):
        """
        Compares the perl reference and maxentpy MES values against each other, using tests/priors_short.tsv only
        for the input sequences.
        """
        self.compare_MSE(compare_type='cross')
