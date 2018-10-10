#!/usr/bin/env python

"""
Description:
    LOVD submissions are differentiated by a combination of the variant (cDNA), submitter, and
    functional analysis technique. Combine the fields to create an id to identify the same submission
    between releases (used for versioning).
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


def main():
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

    # add submission id to each submission (cDNA + submitter + functional analysis technique)
    for line in f_in:
        line = line.replace('"', '')
        parsedLine = line.split('\t')

        cDNA = parsedLine[fieldIdxDict['cDNA']]
        submitter = parsedLine[fieldIdxDict['submitters']]
        functional_analysis_technique = parsedLine[fieldIdxDict['functional_analysis_technique']]

        submissionId = cDNA + submitter + functional_analysis_technique

        ammendedLine = line.split('\r\n')[0] + '\t' + submissionId + '\r\n'
        f_out.write(ammendedLine)


if __name__ == "__main__":
    main()
