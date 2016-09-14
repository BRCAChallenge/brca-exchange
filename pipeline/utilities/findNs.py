#!/usr/bin/env python

import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="built.tsv",
                        help="File with data to check")
    parser.add_argument("--output", default="ns.txt",
                        help="File with ns")

    args = parser.parse_args()
    inputData = csv.DictReader(open(args.input, "r"), delimiter="\t")
    ns = open(args.output, "w")

    # Finds all variants with an N in the column Genomic_Coordinate_hg38
    genomic_coordinates = []
    for row in inputData:
        genomic_coordinates.append(row["Genomic_Coordinate_hg38"])
    sub = "n"
    ns.write("\n".join(s for s in genomic_coordinates if sub.lower() in s.lower()))

if __name__ == "__main__":
    main()
