#!/usr/bin/env python
"""
checks that all bx_ids are accounted for
"""
import argparse
import csv
import logging
from os import listdir
from os.path import isfile, join, abspath
import vcf
import pdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--built",
                        help="built.tsv with all variants")
    parser.add_argument("-o", "--ready_input_dir",
                        help="file directory with all procesed files with bx_ids used to compile built.tsv")
    args = parser.parse_args()

    bx_ids = {}

    files = [f for f in listdir(args.ready_input_dir) if isfile(join(args.ready_input_dir, f)) and "ready" in f or "ENIGMA_combined_with_bx_ids" in f]

    for file in files:
        file_path = abspath(args.ready_input_dir + file)
        if file_path.endswith('.tsv'):
            bx_ids[file] = []
            tsv_file = csv.DictReader(open(file_path, "r"), delimiter='\t')
            for report in tsv_file:
                ids = map(int, report['BX_ID'])
                bx_ids[file] = bx_ids[file] + ids
        else:
            bx_ids[file] = []
            vcf_reader = vcf.Reader(open(file_path, 'r'), strict_whitespace=True)
            try:
                for record in vcf_reader:
                    ids = map(int, record.INFO['BX_ID'])
                    bx_ids[file] = bx_ids[file] + ids
            except ValueError as e:
                print e
    print bx_ids

    csvIn = csv.DictReader(open(args.built, "r"), delimiter='\t')


if __name__ == "__main__":
    main()
