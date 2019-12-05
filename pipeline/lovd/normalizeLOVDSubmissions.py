#!/usr/bin/env python

"""
Description:
    This file constructs a submission id comprised of cdna + submitters + remarks.
    The submission id is later used to combine equivalent submissions.
"""

import argparse
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input LOVD file for ammendment.')
    parser.add_argument('-o', '--output', help='Ouput txt file result.')
    options = parser.parse_args()
    return options

def add_submission_ids(f_in_data_frame):
    for index, row in f_in_data_frame.iterrows():
        # submission id's are used to identify equivalent submissions
        cdna = row['cDNA']
        submitters = row['submitters']
        remarks = row['remarks']
        if isinstance(row['remarks'], basestring):
            remarks = row['remarks']
        else:
            remarks = ''
        f_in_data_frame.loc[index, 'submission_id'] = cdna + submitters + remarks
    return f_in_data_frame


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    f_in_data_frame = pd.read_csv(f_in, sep="\t")
    f_out_data_frame = add_submission_ids(f_in_data_frame)

    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
