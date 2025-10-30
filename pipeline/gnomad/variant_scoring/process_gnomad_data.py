import csv
import gzip
import logging
import shutil
from pathlib import Path

import click

import gnomad.variant_scoring.constants as constants
from common import config as brca_config


def remove_if_exists(path):
    if path.exists():
        shutil.rmtree(path)


def read_and_filter_coverage_v2(input_path, output_path, boundaries):
    """
    Read v2 coverage data from bgzip-compressed TSV and filter line-by-line.

    Args:
        input_path: Path to input .bgz file
        output_path: Path to output CSV file
        boundaries: Dictionary of {chrom: (start, end)} for filtering
    """
    logging.info(f"Reading v2 coverage data from {input_path}")

    with gzip.open(input_path, 'rt') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')

        # Define output columns
        fieldnames = ['chrom', 'pos', 'mean', 'median', 'over_1', 'over_5',
                      'over_10', 'over_15', 'over_20', 'over_25', 'over_30',
                      'over_50', 'over_100']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        rows_written = 0
        for row in reader:
            # Convert chrom to int
            try:
                chrom = int(row['chrom'])
            except (ValueError, KeyError):
                continue

            # Check if within boundaries
            if chrom in boundaries:
                start_bound, end_bound = boundaries[chrom]
                try:
                    pos = int(row['pos'])
                except (ValueError, KeyError):
                    continue

                if start_bound <= pos <= end_bound:
                    # Write filtered row
                    output_row = {field: row.get(field, '') for field in fieldnames}
                    writer.writerow(output_row)
                    rows_written += 1

        logging.info(f"Wrote {rows_written} rows to {output_path}")


def read_and_filter_coverage_v3_v4(input_path, output_path, boundaries):
    """
    Read v3/v4 coverage data and transform to match v2 format, line-by-line.

    Column names and types changed from v2 to v3, so we:
    - Rename 'median_approx' to 'median'
    - Split 'locus' column into 'chrom' and 'pos'
    - Filter out X and Y chromosomes
    - Convert chrom to integer

    Args:
        input_path: Path to input .bgz file
        output_path: Path to output CSV file
        boundaries: Dictionary of {chrom: (start, end)} for filtering
    """
    logging.info(f"Reading v3/v4 coverage data from {input_path}")

    with gzip.open(input_path, 'rt') as infile, open(output_path, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')

        # Define output columns
        fieldnames = ['chrom', 'pos', 'mean', 'median', 'over_1', 'over_5',
                      'over_10', 'over_15', 'over_20', 'over_25', 'over_30',
                      'over_50', 'over_100']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        rows_written = 0
        for row in reader:
            # Split locus column (format: "chr13:12345")
            locus = row.get('locus', '')
            if ':' not in locus:
                continue

            chrom_str, pos_str = locus.split(':', 1)

            # Filter out X and Y chromosomes
            if chrom_str in ['chrX', 'chrY']:
                continue

            # Remove 'chr' prefix and convert to integer
            chrom_str = chrom_str.replace('chr', '')
            try:
                chrom = int(chrom_str)
                pos = int(pos_str)
            except ValueError:
                continue

            # Check if within boundaries
            if chrom in boundaries:
                start_bound, end_bound = boundaries[chrom]
                if start_bound <= pos <= end_bound:
                    # Map median_approx to median
                    output_row = {
                        'chrom': chrom,
                        'pos': pos,
                        'mean': row.get('mean', ''),
                        'median': row.get('median_approx', ''),
                        'over_1': row.get('over_1', ''),
                        'over_5': row.get('over_5', ''),
                        'over_10': row.get('over_10', ''),
                        'over_15': row.get('over_15', ''),
                        'over_20': row.get('over_20', ''),
                        'over_25': row.get('over_25', ''),
                        'over_30': row.get('over_30', ''),
                        'over_50': row.get('over_50', ''),
                        'over_100': row.get('over_100', '')
                    }
                    writer.writerow(output_row)
                    rows_written += 1

        logging.info(f"Wrote {rows_written} rows to {output_path}")


def aggregate_coverage_by_locus(csv_file1, csv_file2, output_path):
    """
    Aggregate two CSV files by locus field and compute combined statistics.

    Processes files line-by-line using CSV module. First pass builds an index,
    second pass merges and writes output.

    Args:
        csv_file1: Path to first CSV file
        csv_file2: Path to second CSV file
        output_path: Path to output CSV file
    """
    logging.info(f"Aggregating coverage data from {csv_file1} and {csv_file2}")

    # First pass: read file1 into a dictionary indexed by locus
    locus_data = {}

    with open(csv_file1, 'r', newline='') as f1:
        reader = csv.DictReader(f1)
        for row in reader:
            locus = f"{row['chrom']}:{row['pos']}"
            try:
                total_dp_1 = float(row.get('mean', '0')) if row.get('mean') else None
                mean_1 = float(row.get('mean', '0')) if row.get('mean') else None
            except (ValueError, TypeError):
                total_dp_1 = None
                mean_1 = None

            locus_data[locus] = {
                'total_DP_1': total_dp_1,
                'mean_1': mean_1,
                'total_DP_2': None,
                'mean_2': None
            }

    # Second pass: read file2 and merge with file1 data
    with open(csv_file2, 'r', newline='') as f2:
        reader = csv.DictReader(f2)
        for row in reader:
            locus = f"{row['chrom']}:{row['pos']}"
            try:
                total_dp_2 = float(row.get('mean', '0')) if row.get('mean') else None
                mean_2 = float(row.get('mean', '0')) if row.get('mean') else None
            except (ValueError, TypeError):
                total_dp_2 = None
                mean_2 = None

            if locus in locus_data:
                locus_data[locus]['total_DP_2'] = total_dp_2
                locus_data[locus]['mean_2'] = mean_2
            else:
                locus_data[locus] = {
                    'total_DP_1': None,
                    'mean_1': None,
                    'total_DP_2': total_dp_2,
                    'mean_2': mean_2
                }

    # Third pass: compute aggregated values and write output
    with open(output_path, 'w', newline='') as outfile:
        fieldnames = ['locus', 'total_DP', 'mean', 'mean_1', 'mean_2', 'total_DP_1', 'total_DP_2']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for locus, data in locus_data.items():
            dp1 = data['total_DP_1']
            dp2 = data['total_DP_2']
            mean1 = data['mean_1']
            mean2 = data['mean_2']

            # Calculate total_DP based on NA conditions
            if dp1 is not None and dp2 is not None:
                total_dp = dp1 + dp2
            elif dp1 is None and dp2 is not None:
                total_dp = dp2
            elif dp1 is not None and dp2 is None:
                total_dp = dp1
            else:
                total_dp = None

            # Calculate mean based on NA conditions
            if mean1 is not None and mean2 is not None:
                if total_dp is not None and total_dp > 0:
                    mean = (mean1 * dp1 + mean2 * dp2) / total_dp
                else:
                    mean = None
            elif mean1 is None and mean2 is not None:
                mean = mean2
            elif mean1 is not None and mean2 is None:
                mean = mean1
            else:
                mean = None

            output_row = {
                'locus': locus,
                'total_DP': total_dp if total_dp is not None else '',
                'mean': mean if mean is not None else '',
                'mean_1': mean1 if mean1 is not None else '',
                'mean_2': mean2 if mean2 is not None else '',
                'total_DP_1': dp1 if dp1 is not None else '',
                'total_DP_2': dp2 if dp2 is not None else ''
            }
            writer.writerow(output_row)

    logging.info(f"Wrote aggregated data to {output_path}")


@click.command()
@click.argument('input_dir', type=click.Path(readable=True))
@click.argument('output_dir', type=click.Path(writable=True))
@click.option('--gene-config-path', type=click.Path(readable=True))
def main(input_dir, output_dir, gene_config_path):
    """
    Process gnomAD coverage data using CSV module for line-by-line processing.

    Reads coverage data from bgzip-compressed TSV files, filters to gene boundaries,
    and outputs CSV files instead of parquet.
    """
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    df_cov_v2_path = output_dir / 'df_cov_v2.csv'
    df_cov_v3_path = output_dir / 'df_cov_v3.csv'
    df_cov_v4_path = output_dir / 'df_cov_v4.csv'

    gene_config = brca_config.load_config(Path(gene_config_path))

    # Build boundaries dictionaries - extract chr column values
    boundaries37 = {}
    boundaries38 = {}

    for idx, row in gene_config.iterrows():
        chrom = row['chr']
        boundaries37[chrom] = (row['start_hg37'], row['end_hg37'])
        boundaries38[chrom] = (row['start_hg38'], row['end_hg38'])

    # Process v2 coverage summaries
    if not df_cov_v2_path.exists():
        logging.info("Processing v2 coverage summaries")
        read_and_filter_coverage_v2(
            input_dir / 'gnomad.exomes.coverage.summary.tsv.bgz',
            df_cov_v2_path,
            boundaries37
        )
    else:
        logging.info(f"Skipping v2 - output already exists: {df_cov_v2_path}")

    # Process v3 coverage summaries
    if not df_cov_v3_path.exists():
        logging.info("Processing v3 coverage summaries")
        read_and_filter_coverage_v3_v4(
            input_dir / 'gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz',
            df_cov_v3_path,
            boundaries38
        )
    else:
        logging.info(f"Skipping v3 - output already exists: {df_cov_v3_path}")

    # Process v4 coverage summaries
    if not df_cov_v4_path.exists():
        logging.info("Processing v4 coverage summaries")
        temp_path = output_dir / 'temp.csv'

        # First read and filter v4 data
        read_and_filter_coverage_v3_v4(
            input_dir / 'gnomad.exomes.v4.0.coverage.summary.tsv.bgz',
            temp_path,
            boundaries38
        )

        # Then aggregate with v3 data
        aggregate_coverage_by_locus(df_cov_v3_path, temp_path, df_cov_v4_path)

        # Clean up temp file
        if temp_path.exists():
            temp_path.unlink()
            logging.info(f"Removed temporary file: {temp_path}")
    else:
        logging.info(f"Skipping v4 - output already exists: {df_cov_v4_path}")

    logging.info("Processing complete!")


if __name__ == "__main__":
    main()
