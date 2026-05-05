#!/usr/bin/env python
"""
Convert a gnomAD coverage CSV file to a parquet file suitable for
popfreq_analysis_prototype_v3.py.

The output parquet contains columns: chrom (int), pos (int), mean (float),
median (float), matching what estimate_coverage() expects.
"""

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert gnomAD coverage CSV to parquet for popfreq_analysis_prototype_v3")
    parser.add_argument("input", help="Input CSV file (e.g. df_cov_v2.csv)")
    parser.add_argument("output", help="Output parquet file (e.g. df_cov_v2.parquet)")
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    df["chrom"] = df["chrom"].astype(int)
    df["pos"] = df["pos"].astype(int)
    df["mean"] = df["mean"].astype(float)
    df["median"] = df["median"].astype(float)
    df = df[["chrom", "pos", "mean", "median"]]
    df.to_parquet(args.output, index=False)
    print(f"Wrote {len(df):,} rows to {args.output}")


if __name__ == "__main__":
    main()
