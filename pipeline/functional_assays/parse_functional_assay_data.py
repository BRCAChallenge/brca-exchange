#!/usr/bin/env python

"""
Description:
    Parses a .tsv source file of functional assays into a list of variants
"""

import argparse
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input functional assay file.')
    parser.add_argument('-o', '--output', help='Ouput tsv file result.')
    options = parser.parse_args()
    return options

def extract_relevant_data_for_processing(f_in_data_frame):

    new_columns = [
        'gene_symbol',
        'reference_sequence',
        'hgvs_nucleotide',
        'log_rna_depletion',
        'functional_enrichment_score']

    new_data_frame = pd.DataFrame(columns=new_columns)

    new_data_frame = new_data_frame.reset_index(drop=True)
    f_in_data_frame = f_in_data_frame.reset_index(drop=True)

    new_data_frame.gene_symbol = f_in_data_frame['gene']
    new_data_frame.reference_sequence = f_in_data_frame['transcript_ID']
    new_data_frame.hgvs_nucleotide = f_in_data_frame['transcript_variant']
    new_data_frame.log_rna_depletion = f_in_data_frame['mean.rna.score']
    new_data_frame.functional_enrichment_score = f_in_data_frame['function.score.mean']

    return new_data_frame


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    # TODO: fixme! quick hack to rem first two lines of extra header categories
    with open(f_in,'r') as f_in:
        with open("/tmp/clean_tmp_functional_assay_data.tsv",'w') as f1:
            # skip header lines
            f_in.next()
            f_in.next()
            for line in f_in:
                f1.write(line)
    f_in = "/tmp/clean_tmp_functional_assay_data.tsv"


    f_in_data_frame = pd.read_csv(f_in, sep="\t")
    f_out_data_frame = extract_relevant_data_for_processing(f_in_data_frame)

    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
