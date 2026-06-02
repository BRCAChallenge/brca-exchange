#!/usr/bin/env python

"""
Description:
    This file constructs a submission id comprised of cdna + submitters + remarks.
    The submission id is later used to combine equivalent submissions.
"""

import argparse
import pandas as pd
import numpy as np

# see https://stackoverflow.com/questions/20625582/how-to-deal-with-settingwithcopywarning-in-pandas
pd.options.mode.chained_assignment = None  # default='warn'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input LOVD file for ammendment.')
    parser.add_argument('-o', '--output', help='Output txt file result.')
    parser.add_argument('-b', '--bracket', help='Bracket variant txt file result.')
    options = parser.parse_args()
    return options

def add_submission_ids(f_in_data_frame):
    for index, row in f_in_data_frame.iterrows():
        # submission id's are used to identify equivalent submissions
        cdna = row['cDNA']
        remarks = row['remarks'] if isinstance(row['remarks'], str) else ''
        submitters = row['submitters'] if isinstance(row['submitters'], str) else ''
        f_in_data_frame.loc[index, 'submission_id'] = cdna + submitters + remarks
    return f_in_data_frame


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output
    bracket_out = options.bracket

    f_in_data_frame = pd.read_csv(f_in, sep="\t")

    # output variants with {} for review
    f_in_data_frame[f_in_data_frame["cDNA"].str.contains("{")==True].to_csv(bracket_out, sep='\t', index=False)

    # remove variants with {} from usable output because downstream hgvs parser can't handle them
    combed_f_in_data_frame = f_in_data_frame[f_in_data_frame["cDNA"].str.contains("{")==False]

    f_out_data_frame = add_submission_ids(combed_f_in_data_frame)

    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
