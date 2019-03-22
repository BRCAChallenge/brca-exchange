#!/usr/bin/env python

"""
Description:
    Parses a .tsv source file of functional assays into a list of variants
"""
import argparse
import pandas as pd
import numpy as np
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input gnomad file.')
    parser.add_argument('-o', '--output', help='Ouput tsv file result.')
    options = parser.parse_args()
    return options

def extract_relevant_data_for_processing(f_in_data_frame):
    relevant_existing_columns = [
        'Chromosome',
        'Position',
        'Reference',
        'Alternate',
        'Source',
        'Allele Count',
        'Allele Number',
        'Allele Frequency',
        'Homozygote Count',
        'Hemizygote Count',
        'Allele Count African',
        'Allele Number African',
        'Homozygote Count African',
        'Hemizygote Count African',
        'Allele Count Latino',
        'Allele Number Latino',
        'Homozygote Count Latino',
        'Hemizygote Count Latino',
        'Allele Count Ashkenazi Jewish',
        'Allele Number Ashkenazi Jewish',
        'Homozygote Count Ashkenazi Jewish',
        'Hemizygote Count Ashkenazi Jewish',
        'Allele Count East Asian',
        'Allele Number East Asian',
        'Homozygote Count East Asian',
        'Hemizygote Count East Asian',
        'Allele Count European (Finnish)',
        'Allele Number European (Finnish)',
        'Homozygote Count European (Finnish)',
        'Hemizygote Count European (Finnish)',
        'Allele Count European (non-Finnish)',
        'Allele Number European (non-Finnish)',
        'Homozygote Count European (non-Finnish)',
        'Hemizygote Count European (non-Finnish)',
        'Allele Count Other',
        'Allele Number Other',
        'Homozygote Count Other',
        'Hemizygote Count Other',
        'Allele Count South Asian',
        'Allele Number South Asian',
        'Homozygote Count South Asian',
        'Hemizygote Count South Asian'
    ]

    column_name_mapping = {
        'Chromosome': 'chr',
        'Position': 'pos_hg19',
        'Reference': 'ref',
        'Alternate': 'alt',
    }

    parsed_data_frame = f_in_data_frame[relevant_existing_columns]
    return parsed_data_frame.rename(columns=column_name_mapping)


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    # NOTE: adjust skiprows as needed to remove extra header lines
    f_in_data_frame = pd.read_csv(f_in)
    f_out_data_frame = extract_relevant_data_for_processing(f_in_data_frame)
    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
