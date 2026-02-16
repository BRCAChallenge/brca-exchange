#!/usr/bin/env python                                                                                                   

"""                                                                                                                      
Description:                                                                                                            
    Takes in a gnomad table and converts it to vcf format.                                                               
"""

import argparse
from pysam import VariantFile
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Parse a gnomAD VCF for the BRCA Exchange pipeline")
    parser.add_argument('-i', '--input', help="Input VCF file, from gnomAD")
    parser.add_argument('-f', '--fourpointone', help="Special functionality for gnomAD version 4.1",
                        default=False)
    parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')
    options = parser.parse_args()
    return options

def main():
    args = parse_args()
    vcf_in = VariantFile(args.input)
    print("ID\tFilters\tLCR")
    for record in vcf_in.fetch():
        all_filters = ""
        delim = ""
        for this_filter in record.filter.keys():
            all_filters = all_filters + delim + this_filter
            delim = ","
        for this_alt in record.alts:
            if args.fourpointone:
                print("%s-%s-%s-%s\t%s" % (re.sub("chr", "", record.chrom),
                                           record.pos, record.ref, this_alt, all_filters))
            else:
                print("%s-%s-%s-%s\t%s\t%d" % (record.chrom, record.pos, record.ref, this_alt,
                                               all_filters, record.info["lcr"]))        
    
if __name__ == "__main__":
    main()
