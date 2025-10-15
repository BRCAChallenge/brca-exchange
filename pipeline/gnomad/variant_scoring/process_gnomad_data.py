import functools
import logging
import operator
import shutil
from pathlib import Path

import click
import pandas as pd

import gnomad.variant_scoring.constants as constants
from common import config as brca_config


def remove_if_exists(path):
    if path.exists():
        shutil.rmtree(path)


def read_raw_coverage_data(input_path):
    """
    Read coverage data from bgzip-compressed TSV file using pandas.

    Pandas can handle bgzip (.bgz) files directly through compression='gzip'.
    """
    return pd.read_csv(str(input_path), sep='\t', compression='gzip')



def boundaries_predicate_variants(boundaries):
    """filter all records within gene boundaries of the variant data"""
    chrom_predicate = functools.reduce(operator.or_,
                                       (((sf.col('contigName') == chrom) & (bound[0] <= sf.col('start')) & (
                                               sf.col('end') <= bound[1])) for (chrom, bound) in boundaries.items()))

    return chrom_predicate


def boundaries_predicate_coverage(boundaries):
    """filter all records within gene boundaries in the variant data"""
    return functools.reduce(operator.or_,
                            (((sf.col('chrom') == chrom) &
                              (bound[0] <= sf.col('pos')) &
                              (sf.col('pos') <= bound[1])) for (chrom, bound) in boundaries.items()))


def read_coverage_data_v2(path):
    """Read v2 coverage data and ensure chrom column is integer type."""
    df_cov = read_raw_coverage_data(str(path))
    df_cov['chrom'] = df_cov['chrom'].astype(int)
    return df_cov


def read_coverage_data_v3_v4(path):
    """
    Read v3/v4 coverage data and transform to match v2 format.

    Column names and types changed from v2 to v3, so we:
    - Rename 'median_approx' to 'median'
    - Split 'locus' column into 'chrom' and 'pos'
    - Filter out X and Y chromosomes
    - Convert chrom to integer
    """
    df = read_raw_coverage_data(str(path))

    # Rename median_approx to median
    df = df.rename(columns={'median_approx': 'median'})

    # Split locus column (format: "chr13:12345") into chrom and pos
    df[['chrom', 'pos']] = df['locus'].str.split(':', expand=True)

    # Filter out X and Y chromosomes
    df = df[~df['chrom'].isin(['chrX', 'chrY'])]

    # Remove 'chr' prefix and convert to integer
    df['chrom'] = df['chrom'].str.replace('chr', '').astype(int)
    df['pos'] = df['pos'].astype(int)

    # Select relevant columns
    columns_to_keep = ['chrom', 'pos', 'mean', 'median', 'over_1', 'over_5',
                       'over_10', 'over_15', 'over_20', 'over_25', 'over_30',
                       'over_50', 'over_100']

    return df[columns_to_keep]



def prepare_coverage_data(coverage_path, coverage_reader, boundaries):
    """
    Prepare coverage data by reading and filtering to gene boundaries.

    Now uses pandas filtering instead of Spark.
    """
    df_cov = coverage_reader(coverage_path)
    # Apply boundaries filter using pandas instead of Spark
    # Build a boolean mask for all boundaries
    mask = pd.Series([False] * len(df_cov), index=df_cov.index)
    for chrom, (start, end) in boundaries.items():
        mask |= ((df_cov['chrom'] == chrom) &
                 (df_cov['pos'] >= start) &
                 (df_cov['pos'] <= end))

    return df_cov[mask]


def aggregate_parquet_by_locus(parquet_file1, parquet_file2):
    """
    Aggregate two parquet files by locus field and sum their total_DP columns.

    Args:
        parquet_file1: Path to first parquet file
        parquet_file2: Path to second parquet file

    Returns:
        pandas.DataFrame with columns: locus, total_DP (sum of both files), mean (weighted average),
        mean_1, mean_2, total_DP_1, total_DP_2
    """
    df1 = pd.read_parquet(parquet_file1)
    df2 = pd.read_parquet(parquet_file2)
    # Select locus, total_DP, and mean columns
    df1_subset = df1[['locus', 'total_DP', 'mean']].rename(columns={'total_DP': 'total_DP_1', 'mean': 'mean_1'})
    df2_subset = df2[['locus', 'total_DP', 'mean']].rename(columns={'total_DP': 'total_DP_2', 'mean': 'mean_2'})
    # Merge on locus
    df_merged = pd.merge(df1_subset, df2_subset, on='locus', how='outer')

    # Calculate total_DP based on NA conditions
    def compute_total_dp(row):
        dp1_is_na = pd.isna(row['total_DP_1'])
        dp2_is_na = pd.isna(row['total_DP_2'])
        if not dp1_is_na and not dp2_is_na:
            # Both not NA: sum them
            return row['total_DP_1'] + row['total_DP_2']
        elif dp1_is_na and not dp2_is_na:
            # Only total_DP_1 is NA: use total_DP_2
            return row['total_DP_2']
        elif not dp1_is_na and dp2_is_na:
            # Only total_DP_2 is NA: use total_DP_1
            return row['total_DP_1']
        else:
            # Both are NA: return NA
            return pd.NA

    df_merged['total_DP'] = df_merged.apply(compute_total_dp, axis=1)

    # Calculate mean based on NA conditions
    def compute_mean(row):
        mean_1_is_na = pd.isna(row['mean_1'])
        mean_2_is_na = pd.isna(row['mean_2'])
        if not mean_1_is_na and not mean_2_is_na:
            # Both not NA: compute weighted average
            if row['total_DP'] > 0:
                return (row['mean_1'] * row['total_DP_1'] + row['mean_2'] * row['total_DP_2']) / row['total_DP']
            else:
                return pd.NA
        elif mean_1_is_na and not mean_2_is_na:
            # Only mean_1 is NA: use mean_2
            return row['mean_2']
        elif not mean_1_is_na and mean_2_is_na:
            # Only mean_2 is NA: use mean_1
            return row['mean_1']
        else:
            # Both are NA: return NA
            return pd.NA

    df_merged['mean'] = df_merged.apply(compute_mean, axis=1)

    # Return locus, total_DP, mean, and intermediate columns
    return df_merged[['locus', 'total_DP', 'mean', 'mean_1', 'mean_2', 'total_DP_1', 'total_DP_2']]


@click.command()
@click.argument('input_dir', type=click.Path(readable=True))
@click.argument('output_dir', type=click.Path(writable=True))
@click.option('--gene-config-path', type=click.Path(readable=True))
def main(input_dir, output_dir, gene_config_path):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    df_cov_v2_path = output_dir + '/df_cov_v2.parquet'
    df_cov_v3_path = output_dir + '/df_cov_v3.parquet'
    df_cov_v4_path = output_dir + '/df_cov_v4.parquet'

    gene_config = brca_config.load_config(Path(gene_config_path))

    boundaries37 = {c: (s, e) for _, (c, s, e) in gene_config[['chr', 'start_hg37', 'end_hg37']].iterrows()}
    boundaries38 = {c: (s, e) for _, (c, s, e) in gene_config[['chr', 'start_hg38', 'end_hg38']].iterrows()}

    if not df_cov_v2_path.exists():
        logging.info("Processing v2 coverage summaries")
        df_cov_v2 = prepare_coverage_data(input_dir + '/gnomad.exomes.coverage.summary.tsv.bgz',
                                          read_coverage_data_v2, boundaries37)
        df_cov_v2.to_parquet(df_cov_v2_path)

    if not df_cov_v3_path.exists():
        logging.info("Processing v3 coverage summaries")
        df_cov_v3 = prepare_coverage_data(input_dir + '/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz',
                                          read_coverage_data_v3_v4, boundaries38)
        df_cov_v3.to_parquet(df_cov_v3_path)

    if not df_cov_v4_path.exists():
        logging.info("Processing v4 coverage summaries")
        df_cov_v4_temp = prepare_coverage_data(input_dir + '/gnomad.exomes.v4.0.coverage.summary.tsv.bgz',
                                               read_coverage_data_v3_v4, boundaries38, spark)
        df_cov_v4_temp.to_parquet(str(output_dir + '/temp.parquet'))
        df_cov_v4 = aggregate_parquet_by_locus(df_cov_v3_path, output_dir + '/temp.parquet')
        df_cov_v4.to_parquet(df_cov_v4_path)


if __name__ == "__main__":
    main()
