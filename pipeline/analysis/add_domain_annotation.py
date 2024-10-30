#!/usr/bin/env python

import argparse
import csv
import re
import sys

DOMAIN_NAME = "Domain"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="build_final.tsv",
                        help="built_final")
    parser.add_argument("-o", "--output", default="built_with_domain.tsv",
                        help="version of input file with new column added")
    parser.add_argument("-b", "--bedfile", default="BRCAclinDomains.bed",
                        help="BED file with the domain definition")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Print debugging info")
    args = parser.parse_args()
    return(args)


def initialize_output_file(input_file, output_filename):
    """
    Create an empty output file with the new columns                        
    """
    new_columns = [DOMAIN_NAME]
    input_header_row = input_file.fieldnames
    if "change_type" in input_header_row:
        idx = input_header_row.index("change_type")
        output_header_row = input_header_row[:idx] + new_columns \
            + input_header_row[idx:]
    else:
        output_header_row = input_header_row + new_columns
    output_file = csv.DictWriter(open(output_filename,"w"),
                                 fieldnames=output_header_row,
                                 delimiter = '\t')
    output_file.writeheader()
    return(output_file)

def read_domain_bed(bedfile):
    domains = {}
    with open(bedfile, 'r') as fp:
        csv_reader = csv.reader(fp, delimiter='\t')
        for row in csv_reader:
            assert(len(row) >= 4)
            interval = {}
            chrom = re.sub("^chr", "", row[0])
            interval["start"] = int(row[1])
            interval["end"] = int(row[2])
            interval["label"] = row[3]
            if not chrom in domains:
                domains[chrom]  = []
            domains[chrom].append(interval)
    return(domains)
            

def find_overlapping_domain(domains, chrom, start, end):
    assert(chrom in domains)
    for interval in domains[chrom]:
        if start <= interval["end"] and end >= (interval["start"]+1):
            return(interval["label"])
    return(None)



    
def main():
    csv.field_size_limit(sys.maxsize)
    args = parse_args()
    domains = read_domain_bed(args.bedfile)
    with open(args.input, 'r') as input_fp:
        input_reader = csv.DictReader(input_fp, delimiter = "\t")
        writer = initialize_output_file(input_reader, args.output)
        for variant in input_reader:
            domain = find_overlapping_domain(domains, variant["Chr"],
                                             int(variant["Hg38_Start"]),
                                             int(variant["Hg38_End"]))
            if domain is None:
                variant[DOMAIN_NAME] = "-"
            else:
                variant[DOMAIN_NAME] = domain
            writer.writerow(variant)
                
if __name__ == "__main__":
    main()
