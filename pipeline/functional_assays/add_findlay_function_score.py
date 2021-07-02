#!/usr/bin/env python

"""
Description:
    Parses a .tsv source file of functional assays into a list of variants
"""
import argparse
import pandas as pd
import numpy as np
import os
import pdb


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--findlay', help='Input findlay functional assay file.')
    parser.add_argument('-e', '--enigma', help='Input enigma functional assay file.')
    parser.add_argument('-o', '--output', help='Ouput tsv file result.')
    options = parser.parse_args()
    return options


def main():
    options = parse_args()
    findlay = options.findlay
    enigma = options.enigma
    f_out = options.output

    # NOTE: adjust skiprows as needed to remove extra header lines
    findlay_data_frame = pd.read_csv(findlay, skiprows=2, sep="\t")
    enigma_data_frame = pd.read_csv(enigma, skiprows=1, sep="\t")

    findlay_data_frame = findlay_data_frame.rename(columns={"transcript_variant": "HGVS Nucleotide Variant"})
    df = enigma_data_frame.merge(findlay_data_frame, on='HGVS Nucleotide Variant', how="left")
    enigma_data_frame['Function Score'] = df['function.score.mean']

    enigma_data_frame.to_csv(options.output, sep='\t', index=False)


if __name__ == "__main__":
    main()
