#!/usr/bin/env python

import argparse
import csv
import glob
import logging
import os
import shutil
import subprocess
import tempfile
from urllib.request import urlretrieve

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

GNOMAD_DOWNLOAD_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/"
    "joint/gnomad.joint.v4.1.sites.chr%s.vcf.bgz"
    )

def parse_args():
    parser = argparse.ArgumentParser(description='Download selected genes to VCF format.')
    parser.add_argument('-g', '--gene_config', type=argparse.FileType('r'),
                        help='Workflow file with data on the selected genes')
    parser.add_argument('-l', '--logfile',
                        default='download_gnomad_fourpointone.log')
    parser.add_argument('-o', '--output', help="output file",
                        default="gnomADv4.hg38.vcf")
    parser.add_argument('-v', '--verbose', action='count', default=False,
                        help='determines logging')
    args = parser.parse_args()
    return args

def process_one_gene(chrom, start_coord, end_coord, output_vcf, logger):
    # Step 1: download the VCF and its .tbi file
    remote_file_vcf = GNOMAD_DOWNLOAD_URL % chrom
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
        subprocess.run(bcftools_cmd, stdout=f_out)
    subprocess.run(["bcftools", "index", "-t", output_vcf])
        
    # Step 3: cleanup.  Remove the big VCF file and its tbi file
    os.remove(local_file_vcf)
    os.remove(local_file_tbi)

    
def main():
    args = parse_args()
    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.CRITICAL
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=args.logfile, filemode="w",
                        level=logging_level)
    gene_config_data = csv.DictReader(args.gene_config)
    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as tmp_fp:
        for row in gene_config_data:
            output_vcf = "%s.gnomADv4.1.joint.vcf.gz" % row["symbol"]
            process_one_gene(row["chr"], row["start_hg38"], row["end_hg38"],
                             output_vcf, logger)
            tmp_fp.write(output_vcf + "\n")
    merge_cmd = ["bcftools", "merge", "--file-list", tmp.name,
                    "-Ov", "-o", args.output]
    subprocess.run(merge_cmd)
    
    #
    # Clean up the individual files and their indices
    for thisfile in glob.glob("*.vcf.gz*"):
        os.remove(thisfile)


    

    
if __name__ == "__main__":
    main()
