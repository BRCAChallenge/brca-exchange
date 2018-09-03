#!/usr/bin/env python

import csv
import unittest
from filterBlacklistedVars import filter_vars_pandas, filter_vars_csv, priorsCols

class TestStringMethods(unittest.TestCase):
    def filterShortPriors(self, method):
        priorsPath = "tests/priors_short.tsv"
        blacklistPath = "tests/blacklisted_vars_test.txt"
        outputPath = "/tmp/priors_short_filtered.tsv"

        # first, produce the file with redacted values
        with open(priorsPath, "r") as variants_fp, \
                open(blacklistPath, "r") as blacklist_fp, \
                open(outputPath, "w") as outfile_fp:
            # gather the blacklist so we can check only those variants later on
            blacklist = set(x.strip() for x in blacklist_fp.readlines())
            blacklist_fp.seek(0)

            # apply whichever method we're testing to filter the variants
            method(variants_fp, blacklist_fp, outfile_fp)

        # check that all the columns in the output file have been redacted
        with open(priorsPath, "r") as originalFP, open(outputPath, "r") as outputFP:
            orig_variants = csv.DictReader(originalFP, delimiter="\t")
            new_variants = csv.DictReader(outputFP, delimiter="\t")

            for old, new in zip(orig_variants, new_variants):
                var_id = "%s:%s" % (old['Reference_Sequence'], old['HGVS_cDNA'])
                for col in old:
                    self.assertEquals(new[col], '-' if var_id in blacklist and col in priorsCols else old[col])

    @unittest.expectedFailure
    def test_filterShortPriorsPandas(self):
        # this test fails because pandas reinterprets floats when it loads a tsv,
        # introducing very small rounding errors
        self.filterShortPriors(method=filter_vars_pandas)

    def test_filterShortPriorsCSV(self):
        self.filterShortPriors(method=filter_vars_csv)
