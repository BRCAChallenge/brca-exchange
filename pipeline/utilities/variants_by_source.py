#!/usr/bin/env python

import argparse
import csv
import collections


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="built.tsv",
                        help="File with data to check")

    num_variants_by_source = {
                              "ClinVar": 0,
                              "ENIGMA": 0,
                              "1000_Genomes": 0,
                              "exLOVD": 0,
                              "LOVD": 0,
                              "ExAC": 0,
                              "BIC": 0,
                              "ESP": 0
                             }

    args = parser.parse_args()
    inputData = csv.DictReader(open(args.input, "r"), delimiter="\t")

    for variant in inputData:
        for source in num_variants_by_source.keys():
            if source in variant["Source"]:
                num_variants_by_source[source] += 1

    print num_variants_by_source

if __name__ == "__main__":
    main()
