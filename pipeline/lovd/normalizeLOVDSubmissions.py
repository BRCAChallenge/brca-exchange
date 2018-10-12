#!/usr/bin/env python

"""
Description:
    LOVD functional analysis technique and result are reported in a colon delimited fashion
    e.g. technique:result. This file splits them into two separate fields, then constructs a
    submission id comprised of cdna + submitters + functional analysis technique. The
    submission id is later used to combine equivalent submissions.
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

def separate_functional_analysis_technique_and_result_and_add_submission_ids(f_in_data_frame):

    # set default empty values
    f_in_data_frame['functional_analysis_technique'] = np.nan
    f_in_data_frame['functional_analysis_result'] = np.nan

    for index, row in f_in_data_frame.iterrows():

        functional_analysis = row['functional_analysis']

        # default used for submission id if no value is present
        functional_analysis_technique = ''

        '''functional_analysis is reported in a colon delimited format,
        technique:result, so we split it into two separate fields'''
        if isinstance(functional_analysis, basestring):
            split_functional_analysis = functional_analysis.split(':', 1)
            functional_analysis_technique = split_functional_analysis[0]
            if len(split_functional_analysis) == 1:
                # no functional_analysis_result value is present
                f_in_data_frame.loc[index, 'functional_analysis_technique'] = functional_analysis_technique
            else:
                f_in_data_frame.loc[index, 'functional_analysis_technique'], f_in_data_frame.loc[index, 'functional_analysis_result'] = split_functional_analysis

        # submission id's are used to identify equivalent submissions
        cdna = row['cDNA']
        submitters = row['submitters']
        f_in_data_frame.loc[index, 'submission_id'] = cdna + submitters + functional_analysis_technique
    return f_in_data_frame


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    f_in_data_frame = pd.read_csv(f_in, sep="\t")
    f_out_data_frame = separate_functional_analysis_technique_and_result_and_add_submission_ids(f_in_data_frame)

    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
