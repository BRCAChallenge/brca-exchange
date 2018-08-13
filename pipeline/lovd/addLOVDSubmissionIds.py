#!/usr/bin/env python

"""
Zack Fischmann
5/8/2018

Description:
    LOVD submissions are differentiated by a combination of the variant (cDNA) and submitter.
    Combine the two fields to create an id to identify the same submission between releases (used
    for versioning).
"""

import argparse
import sys
import os
import urllib
from collections import defaultdict


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input LOVD file for ammendment.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Ouput txt file result.')
    options = parser.parse_args()
    return options


def main(args):
    options = parse_args()
    f_in = options.input
    f_out = options.output

    headerline = f_in.readline()

    # add submission_id field
    ammendedHeaderline = headerline.split('\r\n')[0] + '\t' + '"submission_id"' + '\r\n'
    f_out.write(ammendedHeaderline)

    # get indexes of fields
    parsedHeaderline = headerline.strip().replace(' ', '_').replace('"', '').split('\t')
    fieldIdxDict = defaultdict()
    for index, field in enumerate(parsedHeaderline):
        fieldIdxDict[field] = index

    # add submission id to each submission (cDNA + submitter + functional analysis)
    for line in f_in:
        line = line.replace('"', '')
        parsedLine = line.strip().split('\t')

        cDNA = parsedLine[fieldIdxDict['cDNA']]
        submitter = parsedLine[fieldIdxDict['submitters']]
        functional_analysis = parsedLine[fieldIdxDict['functional_analysis']]

        submissionId = cDNA + submitter + functional_analysis

        ammendedLine = line.split('\r\n')[0] + '\t' + submissionId + '\r\n'
        f_out.write(ammendedLine)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
