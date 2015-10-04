"""this script takes a vcf file, read each row, if the ALT field contains more than 
one item, it will create multiple variant row based on that row
"""

import vcf
import sys
import argparse
from copy import deepcopy

def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    return args


def main():
    args = arg_parse()
    vcf_reader = vcf.Reader(open(args.input, "r"), strict_whitespace=True)
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)

    for record in vcf_reader:
        n = len(record.ALT)
        if n == 1:
            vcf_writer.write_record(record)
        else:
            for i in range(n):
                new_record = deepcopy(record) 
                new_record.ALT = [deepcopy(record.ALT[i])]
                for key in record.INFO.keys():
                    value = deepcopy(record.INFO[key])
                    if type(value) == list and len(value) == n:
                        new_record.INFO[key] = [value[i]]
                vcf_writer.write_record(new_record)

def my_vcf_writer(record):
    



if __name__=="__main__":
    main()
