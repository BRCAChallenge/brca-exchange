#!/usr/bin/env python

import argparse
import csv
import re

def headerLineToColumnIndex(contents):
    """Given a line of header information, return a dictionary indicating
    which column is in which position"""
    columnIndex = dict()
    for ii in range(0, len(contents)):
        columnIndex[contents[ii]] = ii
    return columnIndex

class variantData:

    def __init__(self, pos, ref, alt, significance, hgvs):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.significance = significance
        self.hgvs = hgvs

    def equivalent(self, pos, ref, alt, significance):
        #print self.hgvs, "pos", self.pos, pos, (pos == self.pos)
        #print self.hgvs, "ref", self.ref, ref, (ref == self.ref)
        #print self.hgvs, "alt", self.alt, alt, (alt == self.alt)
        #print self.hgvs, "significance", self.significance, significance, (significance.lower() == self.significance.lower()) 
        if pos != self.pos or ref != self.ref or alt != self.alt \
           or significance.lower() != self.significance.lower():
            #print self.hgvs, "equivalent"
            return False
        else:
            #print self.hgvs, "not equivalent"
            return True
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ucscTabFile")
    parser.add_argument("invitaeVcf")
    args = parser.parse_args()
    
    # Read the UCSC file.  Echo each line, and record the position of each
    # variant by HGVS string.
    variantToPos = dict()
    ucscRowCount = 1
    invitaeSubmissions = dict()
    ucsc = open(args.ucscTabFile)
    for row in ucsc:
        print row.rstrip()
        tokens = row.split("\t")
        if ucscRowCount == 1:
            columnIndex = headerLineToColumnIndex(tokens)
        else:
            variant = tokens[columnIndex["HGVS"]]
            pos = tokens[columnIndex["Pos"]]
            variantToPos[variant] = pos
            submitter = re.sub("(\s)*$", "", tokens[columnIndex["Submitter"]])
            if submitter == "Invitae":
                ref = tokens[columnIndex["Ref"]]
                alt = tokens[columnIndex["Alt"]]
                significance = tokens[columnIndex["ClinicalSignificance"]]
                invitaeSubmissions[variant] = variantData(pos, ref, alt, significance, variant)
        ucscRowCount += 1

    # Read the Invitae coordinates file.  When a variant has been
    # seen before, print it with the previous coordinate.  When it's
    # a new variant, print it with the Invitae coordinate
    invitae = open(args.invitaeVcf)
    invitaeRowCount = 1
    for row in invitae:
        tokens = row.rstrip().split('\t')
        if invitaeRowCount == 1:
            invitaeColumnIndex = headerLineToColumnIndex(tokens)
        else:
            variant = tokens[invitaeColumnIndex["HGVS"]]
            # Check to see if this Invitae variant is already in there...
            alreadySubmitted = False
            significance = tokens[invitaeColumnIndex["interp"]]
            ref = tokens[invitaeColumnIndex["vcf_ref"]]
            alt = tokens[invitaeColumnIndex["vcf_alt"]]
            if invitaeSubmissions.has_key(variant):
                if invitaeSubmissions[variant].equivalent(tokens[invitaeColumnIndex["vcf_pos"]], ref, alt, significance):
                    alreadySubmitted = True
            if not alreadySubmitted:
                if variantToPos.has_key(variant):
                    pos = variantToPos[variant]
                else:
                    pos = tokens[invitaeColumnIndex["vcf_pos"]]
                print '\t'.join((variant, "Invitae", significance, tokens[invitaeColumnIndex["last evaluated"]],
                                 tokens[invitaeColumnIndex["last evaluated"]], "SCV000000000", "NA", 
                                 tokens[invitaeColumnIndex["vcf_chrom"]],
                                 pos, ref, alt, tokens[invitaeColumnIndex["gene"]]))
        invitaeRowCount += 1


if __name__ == "__main__":
    # execute only if run as a script
    main()

                
