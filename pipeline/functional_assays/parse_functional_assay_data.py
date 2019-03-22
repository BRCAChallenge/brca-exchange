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
    parser.add_argument('-i', '--input', help='Input functional assay file.')
    parser.add_argument('-o', '--output', help='Ouput tsv file result.')
    options = parser.parse_args()
    return options

def extract_relevant_data_for_processing(f_in_data_frame):
    relevant_existing_columns = [
        'chromosome',
        'position (hg19)',
        'reference',
        'alt',
        'gene',
        'transcript_ID',
        'transcript_variant',
        'mean.rna.score',
        'function.score.mean'
    ]

    column_name_mapping = {
        'chromosome': 'chr',
        'position (hg19)': 'pos_hg19',
        'reference': 'ref',
        'alt': 'alt',
        'gene': 'gene_symbol',
        'transcript_ID': 'reference_sequence',
        'transcript_variant': 'hgvs_nucleotide',
        'mean.rna.score': 'log_rna_depletion',
        'function.score.mean': 'functional_enrichment_score'
    }

    parsed_data_frame = f_in_data_frame[relevant_existing_columns]
    return parsed_data_frame.rename(columns=column_name_mapping)


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    # NOTE: adjust skiprows as needed to remove extra header lines
    f_in_data_frame = pd.read_csv(f_in, skiprows=2, sep="\t")
    f_out_data_frame = extract_relevant_data_for_processing(f_in_data_frame)
    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
