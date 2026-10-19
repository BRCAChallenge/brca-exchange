#!/Usr/bin/env python

import argparse
import csv
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
    plus genome) or exome-only counts and allele frequencies.  It takes as
    input a gene list, with coordinates, encapsulated in one of the
    gene_config_* files found in the pipeline workflow directory.  For each
    gene, it generates output with the filename '<gene>.gnomADv4.1.<source>.vcf'
"""

GNOMAD_JOINT_DOWNLOAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/"
    "joint/gnomad.joint.v4.1.sites.chr%s.vcf.bgz"
    )
GNOMAD_EXOME_DOWNLOAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/"
    "exomes/gnomad.exomes.v4.1.sites.chr%s.vcf.bgz"
    )

def parse_args():
    parser = argparse.ArgumentParser(description='Download selected genes to VCF format.')
    parser.add_argument('-g', '--gene_config', type=argparse.FileType('r'),
                        help='Workflow file with data on the selected genes')
    parser.add_argument('-l', '--logfile',
                        default='download_gnomad_fourpointone.log')
    parser.add_argument('-o', '--output', help="output file",
                        default="gnomADv4.hg38.vcf")
    parser.add_argument('-c', '--coverage', required=True,
                        help="path to precomputed weighted coverage parquet")
    parser.add_argument('-s', '--source', choices=['joint', 'exome'], default='joint',
                        help="gnomAD v4.1 data source: joint (default) or exome")
    parser.add_argument('-v', '--verbose', action='count', default=False,
                        help='determines logging')
    parser.add_argument('--postprocess-only', action='store_true',
                        help='Skip download; postprocess --input-vcf with coverage and write to --output')
    parser.add_argument('--input-vcf',
                        help='Input VCF for --postprocess-only mode')
    args = parser.parse_args()
    return args

def process_one_gene(chrom, start_coord, end_coord, output_vcf, download_url, logger):
    # Step 1: download the VCF and its .tbi file
    remote_file_vcf = download_url % chrom
    local_file_vcf = chrom + ".vcf"
    logger.info("downloading %s" % remote_file_vcf)
    urlretrieve(remote_file_vcf, local_file_vcf)
    remote_file_tbi = remote_file_vcf + ".tbi"
    local_file_tbi = local_file_vcf + ".tbi"
    urlretrieve(remote_file_tbi, local_file_tbi)

    # Step 2: extract the data for the coordinates needed
    vcf_range = "chr%s%s%s-%s" % (chrom, ":", start_coord, end_coord)
    bcftools_cmd = ["bcftools", "view", "-Oz","-r", vcf_range, local_file_vcf]
    logger.info("Selecting range of %s from %s" % (vcf_range, local_file_vcf))
    with open(output_vcf, "w") as f_out:
        subprocess.run(bcftools_cmd, stdout=f_out, check=True)
    subprocess.run(["bcftools", "index", "-t", output_vcf], check=True)

    # Step 3: cleanup.  Remove the big VCF file and its tbi file
    os.remove(local_file_vcf)
    os.remove(local_file_tbi)


def postprocess(input_vcf, output_vcf, coverage_parquet, logger):
    """Postprocess the gnomAD file by adding variant_id, flags, and coverage.

    Coverage is looked up from the precomputed weighted coverage parquet
    (output of build_weighted_coverage_parquet).  For multi-position variants
    the per-position weighted_mean_coverage values are averaged over the range.
    """
    logger.info("Reading coverage from %s" % coverage_parquet)
    coverage = pd.read_parquet(coverage_parquet)
    cov_col = ('weighted_mean_coverage' if 'weighted_mean_coverage' in coverage.columns
               else 'mean')
    reader = pysam.VariantFile(input_vcf, 'r')
    if 'variant_id' not in reader.header.info:
        reader.header.info.add("variant_id", number=1, type="String",
                               description="gnomAD-style variant ID")
    if 'flags' not in reader.header.info:
        reader.header.info.add("flags", number=1, type="String",
                               description="gnomAD flags, from the FILTER field")
    reader.header.info.add("coverage", number=1, type="String",
                           description="Variant coverage, according to gnomAD")
    for contig in ('13', '17', 'chr13', 'chr17'):
        if contig not in reader.header.contigs:
            reader.header.add_line('##contig=<ID={}>'.format(contig))
    writer = pysam.VariantFile(output_vcf, 'w', header=reader.header)
    for record in reader:
        # Prefer flags already in INFO (source VCF may carry them directly).
        # Fall back to deriving from the FILTER field for VCFs that encode flags there.
        existing_flags = record.info.get('flags')
        has_real_flags = existing_flags is not None and existing_flags != ('-',)
        if not has_real_flags:
            filter_keys = list(record.filter.keys())
            if not filter_keys or filter_keys[0] == "PASS":
                record.info['flags'] = "-"
            else:
                record.info['flags'] = ','.join(filter_keys)
        alt = record.alts[0]
        var_id = "%s-%s-%s-%s" % (re.sub("^chr", "", record.chrom),
                                  str(record.pos), record.ref, alt)
        record.info['variant_id'] = var_id
        chrom_int = int(str(record.chrom).replace('chr', ''))
        variant_start = record.pos
        variant_end = record.pos + len(alt) - len(record.ref)
        region = coverage[
            (coverage['chrom'] == chrom_int) &
            (coverage['pos'] >= variant_start) &
            (coverage['pos'] <= variant_end)
        ]
        cov = region[cov_col].mean() if not region.empty else None
        logger.debug("coverage for %s: %s" % (var_id, cov))
        record.info['coverage'] = str(cov)
        writer.write(record)
    reader.close()
    writer.close()


def main():
    args = parse_args()
    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=args.logfile, filemode="w",
                        level=logging_level)

    if args.postprocess_only:
        if not args.input_vcf:
            raise ValueError("--input-vcf is required with --postprocess-only")
        postprocess(args.input_vcf, args.output, args.coverage, logger)
        return

    # Resolve caller-supplied paths to absolute before changing into the
    # scratch directory below.
    output_path = os.path.abspath(args.output)
    coverage_path = os.path.abspath(args.coverage)

    download_url = (GNOMAD_JOINT_DOWNLOAD_URL if args.source == 'joint'
                    else GNOMAD_EXOME_DOWNLOAD_URL)
    gene_config_data = csv.DictReader(args.gene_config)

    # process_one_gene() and the merge/sort steps below all use bare
    # relative filenames (e.g. "17.vcf", "unsorted.vcf") -- isolate them in
    # a private scratch directory so a concurrent invocation of this script
    # (the Luigi pipeline runs --source joint and --source exome as
    # independent, parallel tasks sharing the same cwd) can't read, write,
    # or delete the same files out from under this one. Placed alongside
    # the output file rather than the system default (often a small root
    # partition) -- a full pre-subsetting chromosome VCF from gnomAD can be
    # many GB, and the pipeline's own output tree is sized for that.
    orig_cwd = os.getcwd()
    with tempfile.TemporaryDirectory(dir=os.path.dirname(output_path)) as scratch_dir:
        os.chdir(scratch_dir)
        try:
            tmp = tempfile.NamedTemporaryFile()
            with open(tmp.name, 'w') as tmp_fp:
                for row in gene_config_data:
                    output_vcf = "%s.gnomADv4.1.%s.vcf.gz" % (row["symbol"], args.source)
                    process_one_gene(row["chr"], row["start_hg38"], row["end_hg38"],
                                     output_vcf, download_url, logger)
                    tmp_fp.write(output_vcf + "\n")
            merge_cmd = ["bcftools", "merge", "--file-list", tmp.name,
                            "-Ov", "-o", "unsorted.vcf"]
            subprocess.run(merge_cmd, check=True)
            sort_cmd = ["bcftools", "sort", "unsorted.vcf",
                        "-Ov", "-o", "sorted.vcf"]
            subprocess.run(sort_cmd, check=True)
            postprocess("sorted.vcf", output_path, coverage_path, logger)
        finally:
            os.chdir(orig_cwd)



if __name__ == "__main__":
    main()
