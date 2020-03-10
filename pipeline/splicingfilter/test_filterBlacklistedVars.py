#!/usr/bin/env python

import os
import csv
import unittest
from .filterBlacklistedVars import filter_vars, priorsCols

class TestBlacklistedVars(unittest.TestCase):
    def test_filterShortPriors(self):
        curPath = os.path.dirname(os.path.realpath(__file__))

        # this test fails because pandas reinterprets floats when it loads a tsv,
        # introducing very small rounding errors
        priorsPath = os.path.join(curPath, "tests", "priors_short.tsv")
        blacklistPath = os.path.join(curPath, "tests", "blacklisted_vars_test.txt")
        outputPath = "/tmp/priors_short_filtered.tsv"

        # first, produce the file with redacted values
        with open(priorsPath, "r") as variants_fp, \
                open(blacklistPath, "r") as blacklist_fp, \
                open(outputPath, "w") as outfile_fp:
            # gather the blacklist so we can check only those variants later on
            blacklist = set(x.strip() for x in blacklist_fp.readlines())
            blacklist_fp.seek(0)

            # apply whichever method we're testing to filter the variants
            filter_vars(variants_fp, blacklist_fp, outfile_fp)

        # check that all the columns in the output file have been redacted
        with open(priorsPath, "r") as originalFP, open(outputPath, "r") as outputFP:
            orig_variants = csv.DictReader(originalFP, delimiter="\t")
            new_variants = csv.DictReader(outputFP, delimiter="\t")

            for old, new in zip(orig_variants, new_variants):
                var_id = "%s:%s" % (old['Reference_Sequence'], old['HGVS_cDNA'])
                for col in old:
                    expected = '-' if var_id in blacklist and col in priorsCols else old[col]

                    try:
                        self.assertAlmostEqual(float(new[col]), float(expected))
                    except ValueError:
                        # it's not a float, use the regular exact equality
                        self.assertEqual(new[col], expected)
