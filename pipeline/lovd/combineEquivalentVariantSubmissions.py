#!/usr/bin/env python

"""
Zack Fischmann
8/14/2018

Description:
    Combine LOVD variants considered to be equivalent submissions.
"""

import argparse
import csv
import logging
from dateutil.parser import parse
import pdb


def mergeRows(oldRow, newRow):
    if oldRow == newRow:
        oldRow['individuals'] = str(int(oldRow['individuals']) + int(newRow['individuals']))
        return oldRow
    else:
        combinedRow = {}
        for key, val in oldRow.iteritems():
            if key == "individuals":
                combinedRow[key] = str(int(val) + int(newRow[key]))
            elif val != newRow[key]:
                if isinstance(val, basestring):
                    oldVal = [val]
                else:
                    oldVal = val
                if isinstance(newRow[key], basestring):
                    newVal = [newRow[key]]
                else:
                    newVal = newRow[key]
                combinedRow[key] = list(set(oldVal + newVal))
            else:
                combinedRow[key] = val
        return combinedRow


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType('r'),
                        help='Input LOVD file.')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        help='Combined output.')
    args = parser.parse_args()

    csvIn = csv.DictReader(args.input, delimiter='\t')
    csvOut = csv.DictWriter(args.output, delimiter='\t',
                            fieldnames=csvIn.fieldnames)
    csvOut.writeheader()

    combinedSubmissions = {}
    for row in csvIn:
        key = row["submission_id"]
        if key not in combinedSubmissions:
            combinedSubmissions[key] = row
        else:
            combinedSubmissions[key] = mergeRows(combinedSubmissions[key], row)
    for submission_id, submission in combinedSubmissions.iteritems():
        # combined values are handled in the pipeline as strings rather than lists
        for key, val in submission.iteritems():
            if not isinstance(val, basestring):
                submission[key] = ','.join(val)
        csvOut.writerow(submission)


if __name__ == "__main__":
    main()

