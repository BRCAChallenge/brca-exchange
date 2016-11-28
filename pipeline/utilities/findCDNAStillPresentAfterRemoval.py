#!/usr/bin/env python

import argparse
import csv

'''
The purpose of this script is to determine which HGVS_cDNA values found in removed.tsv are still
present in the new data. removed.tsv includes variants that are not present in a new release,
but are present in an older release based on the variants' pyhgvs_Genomic_Coordinate_38 value. 
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--newData", default="built.tsv",
                        help="File with data to check")
    parser.add_argument("--removed", default="removed.tsv",
                        help="File with data to check")

    args = parser.parse_args()
    removed = csv.DictReader(open(args.removed, "r"), delimiter="\t")
    newData = csv.DictReader(open(args.newData, "r"), delimiter="\t")

    removedHGVScDNA = []
    removedHGVScDNAStillPresent = []
    removedHGVScDNAStillPresentInHGVScDNACol = []

    # Find all hgvs cdna values from removed variants
    for variant in removed:
        HGVS_cDNA = variant['HGVS_cDNA']
        if HGVS_cDNA is not "-":
            removedHGVScDNA.append(HGVS_cDNA)

    # Find all hgvs_cdna values that are still present in the new data despite being removed
    for variant in newData:
        # Find instances in the actual "HGVS_cDNA" column.
        if variant['HGVS_cDNA'] in removedHGVScDNA:
            removedHGVScDNAStillPresentInHGVScDNACol.append(variant['HGVS_cDNA'])
        # Find instances in any of the columns.
        for prop, value in variant.iteritems():
            for HGVS_cDNA in removedHGVScDNA:
                if HGVS_cDNA in value:
                    removedHGVScDNAStillPresent.append(HGVS_cDNA)

    print "HGVS_cDNA values still in the new data after being removed."
    print removedHGVScDNAStillPresent

    print 'HGVS_cDNA values still in the new data in column "HGVS_cDNA" after being removed.'
    print removedHGVScDNAStillPresentInHGVScDNACol

if __name__ == "__main__":
    main()
