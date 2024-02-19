#!/usr/bin/env python

import argparse
from common import config as brca_config
import pandas as pd
from pathlib import Path

import gnomad.variant_scoring.constants as constants

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_path", default="original_gnomad",
                        help="Input directory")
    parser.add_argument("-o", "--output_path", default="processed_gnomad",
                        help="Output directory")
    parser.add_argument("-c", "--gene_config_path", 
                        help="gene config file for the genes of interest")
    args = parser.parse_args()
    return(args)



def select_gene_subset(df, chrom, start, end):
    df_subset = df.loc[(df['chrom'] == str(chrom)) & (df['pos'] >= start)
                       & (df['pos'] <= end)]
    return(df_subset)


def select_relevant_rows(df, boundaries):
    df_subset = pd.DataFrame
    for this_chrom in boundaries.keys():
        start_pos = boundaries[this_chrom][0]
        end_pos = boundaries[this_chrom][1]
        subset_this_gene = select_gene_subset(df, this_chrom,
                                              start_pos, end_pos)
        if df_subset.empty:
            df_subset = subset_this_gene
        else:
            df_subset = df_subset.append(subset_this_gene)
    return(df_subset)


def save_relevant_rows(input_filename, boundaries, output_filename):
    df = pd.read_csv(input_filename, sep='\t', dtype={"chrom":str})
    relevant_rows = select_relevant_rows(df, boundaries)
    relevant_rows.to_parquet(output_filename)


def main():
    args = parse_args()
    gene_config = brca_config.load_config(Path(args.gene_config_path))
    boundaries37 = {c: (s, e) for _, (c, s, e) in gene_config[['chr', 'start_hg37', 'end_hg37']].iterrows()}
    boundaries38 = {c: (s, e) for _, (c, s, e) in gene_config[['chr', 'start_hg38', 'end_hg38']].iterrows()}
    #
    # Process the hg37 exomes
    #exomes_37_in = args.input_path + "/gnomad.v2.exomes.coverage.summary.tsv"
    #exomes_37_out = args.output_path + "df_cov_v2.exome.parquet"
    #save_relevant_rows(exomes_37_in, boundaries37, exomes_37_out) 
    #
    # Process the hg37 genomes
    genomes_37_in = args.input_path + "/gnomad.v2.genomes.coverage.summary.tsv"
    genomes_37_out = args.output_path + "df_cov_v2.genomes.parquet"
    save_relevant_rows(genomes_37_in, boundaries37, genomes_37_out) 
    #
    # Process the hg38 genomes
    genomes_38_in = args.input_path + "/gnomad.v3.genomes.coverage.summary.tsv"
    genomes_38_out = args.output_path + "df_cov_v3.genomes.parquet"
    save_relevant_rows(genomes_38_in, boundaries38, genomes_38_out) 


if __name__ == "__main__":
    main()
