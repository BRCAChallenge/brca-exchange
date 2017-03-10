#!/usr/bin/env python
"""
adds bx_id for all enigma reports
"""
import argparse
import csv
import logging


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="input file with merged ENIGMA data")
    parser.add_argument("-o", "--output",
                        help="Output file with corrected ENIGMA data")
    args = parser.parse_args()

    csvIn = csv.DictReader(open(args.input, "r"), delimiter='\t')
    fieldnames = csvIn.fieldnames
    fieldnames.append('BX_ID')
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter='\t',
                            fieldnames=fieldnames)
    csvOut.writerow(dict((fn, fn) for fn in csvIn.fieldnames))
    count = 1
    for row in csvIn:
        row['BX_ID'] = count
        count += 1
        csvOut.writerow(row)


if __name__ == "__main__":
    main()
