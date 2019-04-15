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
        'genome_AFR_ac',
        'genome_AFR_ac_hemi',
        'genome_AFR_ac_hom',
        'genome_AFR_an',
        'genome_AMR_ac',
        'genome_AMR_ac_hemi',
        'genome_AMR_ac_hom',
        'genome_AMR_an',
        'genome_ASJ_ac',
        'genome_ASJ_ac_hemi',
        'genome_ASJ_ac_hom',
        'genome_ASJ_an',
        'genome_EAS_ac',
        'genome_EAS_ac_hemi',
        'genome_EAS_ac_hom',
        'genome_EAS_an',
        'genome_FIN_ac',
        'genome_FIN_ac_hemi',
        'genome_FIN_ac_hom',
        'genome_FIN_an',
        'genome_NFE_ac',
        'genome_NFE_ac_hemi',
        'genome_NFE_ac_hom',
        'genome_NFE_an',
        'genome_OTH_ac',
        'genome_OTH_ac_hemi',
        'genome_OTH_ac_hom',
        'genome_OTH_an',
        'genome_SAS_ac',
        'genome_SAS_ac_hemi',
        'genome_SAS_ac_hom',
        'genome_SAS_an',
        'genome_ac',
        'genome_an',
        'genome_af',
        'exome_AFR_ac',
        'exome_AFR_ac_hemi',
        'exome_AFR_ac_hom',
        'exome_AFR_an',
        'exome_AMR_ac',
        'exome_AMR_ac_hemi',
        'exome_AMR_ac_hom',
        'exome_AMR_an',
        'exome_ASJ_ac',
        'exome_ASJ_ac_hemi',
        'exome_ASJ_ac_hom',
        'exome_ASJ_an',
        'exome_EAS_ac',
        'exome_EAS_ac_hemi',
        'exome_EAS_ac_hom',
        'exome_EAS_an',
        'exome_FIN_ac',
        'exome_FIN_ac_hemi',
        'exome_FIN_ac_hom',
        'exome_FIN_an',
        'exome_NFE_ac',
        'exome_NFE_ac_hemi',
        'exome_NFE_ac_hom',
        'exome_NFE_an',
        'exome_OTH_ac',
        'exome_OTH_ac_hemi',
        'exome_OTH_ac_hom',
        'exome_OTH_an',
        'exome_SAS_ac',
        'exome_SAS_ac_hemi',
        'exome_SAS_ac_hom',
        'exome_SAS_an',
        'exome_an',
        'exome_ac',
        'exome_af',
        'alt',
        'chrom',
        'hgvsc',
        'hgvsp',
        'pos',
        'ref',
        'variantId',
        'consequence'
    ]

    column_name_mapping = {
        'chrom': 'chr',
        'pos': 'pos_hg19',
    }

    parsed_data_frame = f_in_data_frame[relevant_existing_columns]
    return parsed_data_frame.rename(columns=column_name_mapping)


def compile_allele_values(df):
    populations = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'EAS', 'NFE', 'OTH', 'SAS']
    # try:
    df['ac'] = add_values(df['genome_ac'], df['exome_ac'])
    df['an'] = add_values(df['genome_an'], df['exome_an'])
    df['af'] = calculate_frequency(df['ac'], df['an'])
    for population in populations:
        df[population + '_ac'] = add_values(df['genome_' + population + '_ac'], df['exome_' + population + '_ac'])
        df[population + '_ac_hom'] = add_values(df['genome_' + population + '_ac_hom'], df['exome_' + population + '_ac_hom'])
        df[population + '_ac_hemi'] = add_values(df['genome_' + population + '_ac_hemi'], df['exome_' + population + '_ac_hemi'])
        df[population + '_an'] = add_values(df['genome_' + population + '_an'], df['exome_' + population + '_an'])
    df[population + '_af'] = calculate_frequency(df[population + '_ac'], df[population + '_an'])
    return df


def add_values(field_one, field_two):
    return pd.to_numeric(field_one, errors='coerce').add(pd.to_numeric(field_two, errors='coerce'))


def calculate_frequency(ac, an):
    return pd.to_numeric(ac, errors='coerce').divide(pd.to_numeric(an, errors='coerce'))


def main():
    options = parse_args()
    f_in = options.input
    f_out = options.output

    f_in_data_frame = pd.read_csv(f_in, sep='\t')
    relevant_fields_data_frame = extract_relevant_data_for_processing(f_in_data_frame)
    f_out_data_frame = compile_allele_values(relevant_fields_data_frame)
    f_out_data_frame.to_csv(f_out, sep='\t', index=False)


if __name__ == "__main__":
    main()
