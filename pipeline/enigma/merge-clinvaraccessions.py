#!/usr/bin/env python

import argparse
import csv
import re
import os
import sys


def main():
    '''
    This script merges clinvar accession numbers from a .txt file provided by clinvar into enigma data provided by
    enigma. Note that the Enigma file is usually provided as a .xlsx file and must be converted using the steps
    outlined in ./README.md.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--enigma_input",
                        help="Enigma variants")
    parser.add_argument("-c", "--clinvar_input",
                        help="Clinvar accessions")
    parser.add_argument("-o", "--output_file",
                        help="Merged accessions file")

    args = parser.parse_args()

    if '.tsv' not in args.output_file:
        sys.exit("Output file must be of type .tsv")

    final_output_tsv = args.output_file
    temp_csv = args.output_file.replace('.tsv', '.csv')

    #
    # The ClinVar accessions file will begin with a header, with pound signs in the first
    # column.  Most annoyingly, the header row also contains pound signs in the first
    # column, so we can't simply ignore lines with pound signs.  The header row is the
    # first line with tab-delimited content.  So, read down to the first line with
    # more than a single token, and interpret that as the header, and then continue
    # to interpret the file as a normal csv.DictReader.
    #
    fp = open(args.clinvar_input, "rb")
    past_header = False
    while not past_header:
        row = fp.next()
        tokens = row.split('\t')
        if re.search("Your_variant_id\t", row):
            header = tokens
            past_header = True
    clinvar_accessions = csv.DictReader(fp, fieldnames=header, delimiter='\t')
    enigma_variants = csv.DictReader(open(args.enigma_input, "rb"), delimiter="\t")

    with open(temp_csv, 'wb') as csvfile:
        fieldnames = enigma_variants.fieldnames
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

        for clinvar_variant in clinvar_accessions:
            hgvs = clinvar_variant['Your_variant_description'].split(":")[1]
            accession = clinvar_variant['SCV']
            print hgvs, accession

            enigma_variants = csv.DictReader(open(args.enigma_input, "rb"), delimiter="\t")
            for enigma_variant in enigma_variants:
                if enigma_variant['HGVS'] == hgvs:
                    enigma_variant['ClinVarAccession'] = accession
                    writer.writerow(enigma_variant)

    # convert csv output to tsv for final output
    csv.writer(file(final_output_tsv, 'w+'), delimiter="\t").writerows(csv.reader(open(temp_csv)))
    os.remove(temp_csv)

if __name__ == "__main__":
    main()
