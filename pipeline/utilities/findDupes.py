#!/usr/bin/env python

import argparse
import csv
import collections


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="built.tsv",
                        help="File with data to check")
    parser.add_argument("--output", default="dupes.txt",
                        help="File with duplicates")

    args = parser.parse_args()
    inputData = csv.DictReader(open(args.input, "r"), delimiter="\t")
    duplicates = open(args.output, "w")

    # Finds all duplicate variants in a given tsv file containing rows with the column Genomic_Coordinate_hg38
    hgvs_strings = []
    for row in inputData:
        hgvs_strings.append(row["Genomic_Coordinate_hg38"])
    c = collections.Counter(hgvs_strings)
    for hgvs_string, count in c.iteritems():
        if count > 1:
            duplicates.write(hgvs_string + ": " + str(count) + "\n")

if __name__ == "__main__":
    main()
