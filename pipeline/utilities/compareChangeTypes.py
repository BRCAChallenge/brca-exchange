#!/usr/bin/env python

import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--old", default="built_with_change_types.tsv",
                        help="First version of file")
    parser.add_argument("--new", default="built_with_change_types.tsv",
                        help="Second version of file")

    num_variants_by_change_type = {
                                    "new": 0,
                                    "changed_classification": 0,
                                    "changed_information": 0,
                                    "added_information": 0,
                                  }

    args = parser.parse_args()
    oldData = csv.DictReader(open(args.old, "r"), delimiter="\t")
    newData = csv.DictReader(open(args.new, "r"), delimiter="\t")

    variants = {}
    variants_with_incorrect_change_types = {}

    for variant in newData:
        genomic_coordinate = variant["pyhgvs_Genomic_Coordinate_38"]
        change_type = variant["change_type"]
        variants[genomic_coordinate] = change_type

    for variant in oldData:
        genomic_coordinate = variant["pyhgvs_Genomic_Coordinate_38"]
        change_type = variant["change_type"]
        if variants[genomic_coordinate] != change_type:
            variants_with_incorrect_change_types[genomic_coordinate] = [change_type, variants[genomic_coordinate]]

    print variants_with_incorrect_change_types
    print len(variants_with_incorrect_change_types.keys())

if __name__ == "__main__":
    main()
