#!/usr/bin/env python
"""
Find_missing_vcf_records: given VCF files A and B, report the variants in B
which are not in A, in a literal string match of CHR, POS, REF and ALT.  
Output a VCF file with the variants from A that were not found in B.
The output file has the header and the INFO fields of A 
"""

import argparse
import vcf

def variant_id_string(record):
    variant = "%s:%s:%s:%s" % (record.CHROM, str(record.POS),
                               record.REF, record.ALT)
    return(variant)
    

def read_variants_from_b(vcf_b_filename):
    variants_in_b = {}
    reader = vcf.Reader(open(vcf_b_filename, "r"))
    for record in reader:
        variants_in_b[variant_id_string(record)] = 1
    return(variants_in_b)
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--vcf_a",
                        help="VCF with variants that might not be in B")
    parser.add_argument("-b", "--vcf_b",
                        help="VCF that might be missing variants found in A")
    parser.add_argument("-o", "--vcf_out",
                        help="Output VCF file")
    args = parser.parse_args()
    record_count = 0
    vcf_b_variants = read_variants_from_b(args.vcf_b)
    reader = vcf.Reader(open(args.vcf_a, "r")) 
    writer = vcf.Writer(open(args.vcf_out, "w"), reader)
    for record in reader:
        variant = variant_id_string(record)
        if not variant in vcf_b_variants:
            writer.write_record(record)
            record_count += 1
    print(record_count, "records written")


if __name__ == "__main__":
    main()


