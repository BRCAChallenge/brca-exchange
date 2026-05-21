#!/Usr/bin/env python

import argparse
from common import config
import glob
import logging
import os
import pandas as pd
import re
import subprocess
import tempfile
from urllib.request import urlretrieve
import pysam

"""
Description:
    Downloads the gnomAD v4.1 data for the selected genes and translates
    them to VCF.

    This script generates one VCF per selected gene, and the joint (exome 
    plus genome) counts and allele frequencies.  It takes as input a gene
    list, with coordinates, encapsulated in one of the gene_config_* files
    found in the pipeline workflow directory.  For each gene, it generates
    output with the filename '<gene>.gnomADv4.1.vcf'
"""

GNOMAD_VCF_DOWNLOAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/"
    "joint/gnomad.joint.v4.1.sites.chr%s.vcf.bgz"
    )
_PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GNOMAD_GENOMES_DOWNLOAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/"
    "coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"
    )
GNOMAD_EXOMES_DOWNLOAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/"
    "coverage/exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz"
    )
GNOMAD_GENOME_COVERAGE_TSV = os.path.join(_PIPELINE_DIR, "coverage.genome.tsv.bgz")
GNOMAD_EXOME_COVERAGE_TSV = os.path.join(_PIPELINE_DIR, "coverage.exome.tsv.bgz")
GNOMAD_GENOME_COVERAGE_PARQUET = os.path.join(_PIPELINE_DIR, "coverage.genome.parquet")
GNOMAD_EXOME_COVERAGE_PARQUET = os.path.join(_PIPELINE_DIR, "coverage.exome.parquet")
GNOMAD_COMBINED_COVERAGE_PARQUET = os.path.join(_PIPELINE_DIR, "coverage.combined.parquet")


def parse_args():
    parser = argparse.ArgumentParser(description='Download selected genes to VCF format.')
    parser.add_argument('-g', '--gene_config', type=argparse.FileType('r'),
                        help='Workflow file with data on the selected genes')
    parser.add_argument('-l', '--logfile',
                        default='download_gnomad_fourpointone.log')
    parser.add_argument('-o', '--output', help="output file",
                        default="gnomADv4.hg38.vcf")
    parser.add_argument('-c', '--coverage_output', help="output parquet file for weighted coverage",
                        default=GNOMAD_COMBINED_COVERAGE_PARQUET)
    parser.add_argument('-v', '--verbose', action='count', default=False,
                        help='determines logging')
    args = parser.parse_args()
    return args

def process_one_gene(chrom, start_coord, end_coord, output_vcf, logger):
    # Step 1: download the VCF and its .tbi file
    remote_file_vcf = GNOMAD_VCF_DOWNLOAD_URL % chrom
    local_file_vcf = str(chrom) + ".vcf"
    logger.info("downloading VCF file %s" % remote_file_vcf)
    urlretrieve(remote_file_vcf, local_file_vcf)
    remote_file_tbi = remote_file_vcf + ".tbi"
    local_file_tbi = local_file_vcf + ".tbi"
    urlretrieve(remote_file_tbi, local_file_tbi)

    # Step 2: extract the data for the coordinates needed
    vcf_range = "chr%s%s%s-%s" % (chrom, ":", start_coord, end_coord)
    bcftools_cmd = ["bcftools", "view", "-Oz","-r", vcf_range, local_file_vcf]
    logger.info("Selecting range of %s from %s" % (vcf_range, local_file_vcf))
    with open(output_vcf, "w") as f_out:
        subprocess.run(bcftools_cmd, stdout=f_out)
    subprocess.run(["bcftools", "index", "-t", output_vcf])
        
    # Step 3: cleanup.  Remove the big VCF file and its tbi file
    os.remove(local_file_vcf)
    os.remove(local_file_tbi)


def read_coverage(tsv_path, gene_boundaries):
    """
    Read mean coverage and total_DP for the given gene boundaries from a
    local bgzipped gnomAD coverage TSV file, returning a pandas DataFrame
    with columns 'chrom' (int, no 'chr' prefix), 'pos' (int), 'mean' (float),
    and 'total_DP' (float).

    The input data is expected to have a 'locus' column formatted as 'chrN:POS',
    a 'mean' column for mean read depth, and a 'total_DP' column for total
    depth of coverage.
    """
    chrom_ranges = {
        'chr' + str(chrom): (int(str(chrom).replace('chr', '')), start, end)
        for chrom, (start, end) in gene_boundaries.items()
    }
    rows = []
    for chunk in pd.read_csv(tsv_path, sep='\t', compression='gzip',
                              usecols=['locus', 'mean', 'total_DP'],
                              chunksize=100_000):
        chunk[['chrom_str', 'pos_str']] = chunk['locus'].str.split(':', n=1, expand=True)
        chunk['pos'] = chunk['pos_str'].astype(int)
        for chrom_str, (chrom_int, start, end) in chrom_ranges.items():
            mask = ((chunk['chrom_str'] == chrom_str) &
                    (chunk['pos'] >= start) &
                    (chunk['pos'] <= end))
            subset = chunk.loc[mask, ['pos', 'mean', 'total_DP']].copy()
            if not subset.empty:
                subset['chrom'] = chrom_int
                rows.append(subset[['chrom', 'pos', 'mean', 'total_DP']])
    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame(columns=['chrom', 'pos', 'mean', 'total_DP'])


def build_coverage_parquet(tsv_path, gene_boundaries, parquet_path, logger):
    """
    Read coverage data for the given gene boundaries from one gnomAD coverage
    TSV file, write it to a parquet file, and discard the in-memory DataFrame.

    Args:
        tsv_path:       Path to the bgzipped gnomAD coverage TSV.
        gene_boundaries: Dict mapping chrom -> (start, end).
        parquet_path:   Destination path for the parquet file.
        logger:         Logger instance.

    Returns:
        parquet_path
    """
    logger.info("Reading coverage from %s" % tsv_path)
    df = read_coverage(tsv_path, gene_boundaries)
    logger.info("Saving coverage to %s" % parquet_path)
    df.to_parquet(parquet_path, index=False)
    del df
    return parquet_path


def prepare_coverage_data(logger, gene_boundaries):
    """
    Download the gnomAD genome and exome coverage TSV files if needed, then
    convert each to a parquet file.  Each DataFrame is written to disk and
    removed from memory before the next file is read, keeping peak memory use
    to one coverage dataset at a time.

    Returns:
        Tuple of (genome_parquet_path, exome_parquet_path).
    """
    if not os.path.exists(GNOMAD_GENOME_COVERAGE_TSV):
        logger.info("Downloading genome coverage from %s" % GNOMAD_GENOMES_DOWNLOAD_URL)
        urlretrieve(GNOMAD_GENOMES_DOWNLOAD_URL, GNOMAD_GENOME_COVERAGE_TSV)
    if not os.path.exists(GNOMAD_EXOME_COVERAGE_TSV):
        logger.info("Downloading exome coverage from %s" % GNOMAD_EXOMES_DOWNLOAD_URL)
        urlretrieve(GNOMAD_EXOMES_DOWNLOAD_URL, GNOMAD_EXOME_COVERAGE_TSV)
    genome_parquet = build_coverage_parquet(
        GNOMAD_GENOME_COVERAGE_TSV, gene_boundaries, GNOMAD_GENOME_COVERAGE_PARQUET, logger)
    exome_parquet = build_coverage_parquet(
        GNOMAD_EXOME_COVERAGE_TSV, gene_boundaries, GNOMAD_EXOME_COVERAGE_PARQUET, logger)
    return genome_parquet, exome_parquet

    
def build_weighted_coverage_parquet(genome_parquet, exome_parquet, output_path, logger):
    """
    Compute the total_DP-weighted average coverage for every position present
    in either the genome or exome coverage parquets, and write the result to
    a new parquet file.

    Positions present in only one dataset are weighted entirely by that
    dataset's total_DP.  Positions where both total_DP values are zero are
    written with a NaN weighted_mean_coverage.

    Output columns: chrom (int), pos (int), weighted_mean_coverage (float).

    Args:
        genome_parquet: Path to the genome coverage parquet.
        exome_parquet:  Path to the exome coverage parquet.
        output_path:    Destination path for the combined parquet.
        logger:         Logger instance.

    Returns:
        output_path
    """
    logger.info("Reading genome coverage from %s" % genome_parquet)
    genome_df = pd.read_parquet(genome_parquet, columns=['chrom', 'pos', 'mean', 'total_DP'])
    logger.info("Reading exome coverage from %s" % exome_parquet)
    exome_df = pd.read_parquet(exome_parquet, columns=['chrom', 'pos', 'mean', 'total_DP'])
    merged = genome_df.merge(exome_df, on=['chrom', 'pos'], how='outer',
                             suffixes=('_genome', '_exome'))
    for col in ('mean_genome', 'mean_exome', 'total_DP_genome', 'total_DP_exome'):
        merged[col] = merged[col].fillna(0.0)
    total_dp = merged['total_DP_genome'] + merged['total_DP_exome']
    weighted = ((merged['mean_genome'] * merged['total_DP_genome'] +
                 merged['mean_exome'] * merged['total_DP_exome']) / total_dp)
    result = merged[['chrom', 'pos']].copy()
    result['weighted_mean_coverage'] = weighted
    logger.info("Saving combined coverage to %s" % output_path)
    result.to_parquet(output_path, index=False)
    return output_path



    
def main():
    args = parse_args()
    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=args.logfile, filemode="w",
                        level=logging_level)
    
    gene_config_data = config.load_config(args.gene_config)
    boundaries38 = {c: (s, e) for _, (c, s, e)
                    in gene_config_data[['chr', 'start_hg38', 'end_hg38']].iterrows()}
    (genome_parquet,
     exome_parquet) = prepare_coverage_data(logger, boundaries38)
    build_weighted_coverage_parquet(genome_parquet, exome_parquet,
                                    args.coverage_output, logger)
    for f in (GNOMAD_GENOME_COVERAGE_TSV, GNOMAD_EXOME_COVERAGE_TSV,
              genome_parquet, exome_parquet):
        os.remove(f)


if __name__ == "__main__":
    main()
