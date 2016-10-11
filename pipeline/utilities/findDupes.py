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
    parser.add_argument("--column", default="pyhgvs_Genomic_Coordinate_38",
                        help="Column to search for duplicates")

    args = parser.parse_args()
    inputData = csv.DictReader(open(args.input, "r"), delimiter="\t")
    duplicates = open(args.output, "w")

    # Finds all duplicate variants in a given tsv file containing rows with the given column
    instances = []
    for row in inputData:
        instances.append(row[args.column])
    c = collections.Counter(instances)
    for instance, count in c.iteritems():
        if count > 1:
            duplicates.write(instance + ": " + str(count) + "\n")

if __name__ == "__main__":
    main()
