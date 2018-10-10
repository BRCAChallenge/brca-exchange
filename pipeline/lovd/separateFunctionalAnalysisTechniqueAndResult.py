#!/usr/bin/env python

"""
Description:
    LOVD functional analysis technique and result are reported in a colon delimited fashion
    e.g. technique:result. This file splits them into two separate fields.
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

    # add functional analysis fields
    ammendedHeaderline = headerline.split('\r\n')[0] + '\t' + '"functional_analysis_technique"' + '\t' + '"functional_analysis_result"' + '\r\n'
    f_out.write(ammendedHeaderline)

    # get indexes of fields
    parsedHeaderline = headerline.strip().replace(' ', '_').replace('"', '').split('\t')
    fieldIdxDict = defaultdict()
    for index, field in enumerate(parsedHeaderline):
        fieldIdxDict[field] = index

    # fill out data for fields
    for line in f_in:
        line = line.replace('"', '')
        parsedLine = line.strip().split('\t')


        functional_analysis = parsedLine[fieldIdxDict['functional_analysis']]
        functional_analysis_technique = ''
        functional_analysis_result = ''

        split_functional_analysis = functional_analysis.split(':', 1)
        if len(split_functional_analysis) == 1:
            functional_analysis_technique = split_functional_analysis[0]
        else:
            functional_analysis_technique, functional_analysis_result = split_functional_analysis


        ammendedLine = line.split('\r\n')[0] + '\t' + functional_analysis_technique + '\t' + functional_analysis_result + '\r\n'
        f_out.write(ammendedLine)


if __name__ == "__main__":
    main()
