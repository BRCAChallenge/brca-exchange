#!/usr/bin/env python
"""
From ESP, extract the EA (european ancestry) or AA (african american ancestry)
MAFs for the BRCA1 and BRCA2 regions.  produce data in the following format:

chrom_pos_ref_alt Freq

e.g. 13_32890339_G_C0.000199681
"""
import argparse
import vcf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputVcf")
    parser.add_argument("-s", "--start")
    parser.add_argument("-e", "--end")
    parser.add_argument("-a", "--ancestry")
    args = parser.parse_args()

    start = int(args.start)
    end = int(args.end)
    reader = vcf.Reader(open(args.inputVcf, 'r'))
    for record in reader:
        if int(record.POS) >= start and int(record.POS) <= end:
            if args.ancestry == "EA":
                maf = record.INFO["MAF"][0]
            elif args.ancestry == "AA":
                maf = record.INFO["MAF"][1]
            for alt in record.ALT:
                print "%s_%s_%s_%s %s" % (record.CHROM, record.POS, record.REF,
                                          alt, maf)


if __name__ == "__main__":
    # execute only if run as a script
    main()

            
