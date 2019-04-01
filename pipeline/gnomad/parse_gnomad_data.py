#!/usr/bin/env python

"""
Description:
    Parses a .tsv source file of gnomAD data
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
        'AFR_ac',
        'AFR_ac_hemi',
        'AFR_ac_hom',
        'AFR_an',
        'AMR_ac',
        'AMR_ac_hemi',
        'AMR_ac_hom',
        'AMR_an',
        'ASJ_ac',
        'ASJ_ac_hemi',
        'ASJ_ac_hom',
        'ASJ_an',
        'EAS_ac',
        'EAS_ac_hemi',
        'EAS_ac_hom',
        'EAS_an',
        'FIN_ac',
        'FIN_ac_hemi',
        'FIN_ac_hom',
        'FIN_an',
        'NFE_ac',
        'NFE_ac_hemi',
        'NFE_ac_hom',
        'NFE_an',
        'OTH_ac',
        'OTH_ac_hemi',
        'OTH_ac_hom',
        'OTH_an',
        'SAS_ac',
        'SAS_ac_hemi',
        'SAS_ac_hom',
        'SAS_an',
        'ac',
        'ac_hemi',
        'ac_hom',
        'af',
        'alt',
        'an',
        'chrom',
        'datasets',
        'hgvsc',
        'hgvsp',
        'pos',
        'ref',
        'variantId',
    ]

    column_name_mapping = {
        'chrom': 'chr',
        'pos': 'pos_hg19',
    }

    parsed_data_frame = f_in_data_frame[relevant_existing_columns]
    return parsed_data_frame.rename(columns=column_name_mapping)


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    f_in_data_frame = pd.read_csv(f_in, sep='\t')
    f_out_data_frame = extract_relevant_data_for_processing(f_in_data_frame)
    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
