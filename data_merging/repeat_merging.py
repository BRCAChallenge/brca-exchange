"""
this script takes a vcf file, collapses representative variant rows and write out
to a new vcf file (without header)
"""

import vcf
import argparse
import copy
import sys


def arg_parse():
    parser = argparse.ArgumentParser() 
    parser.add_argument("-i", "--input") 
    parser.add_argument("-o", "--output") 
    args = parser.parse_args() 
    return args 

def main():
    print "warning: this script doesn't write out headers"
    args = arg_parse()
    vcf_reader = vcf.Reader(open(args.input, "r"), strict_whitespace=True)
    variant_dict = {}
    num_repeats = 0
    for record in vcf_reader:
        genome_coor = "chr{0}:{1}:{2}>{3}".format(
            record.CHROM, str(record.POS), record.REF, record.ALT[0])
        if genome_coor not in variant_dict.keys():
            variant_dict[genome_coor] = copy.deepcopy(record)
        else:
            num_repeats += 1
            for key in record.INFO:        
                if key not in variant_dict[genome_coor].INFO.keys():
                    variant_dict[genome_coor].INFO[key] = copy.deepcopy(record.INFO[key])
                else:
                    new_value = copy.deepcopy(record.INFO[key])
                    old_value = copy.deepcopy(variant_dict[genome_coor].INFO[key])

                    if type(new_value) != list:
                        new_value = [new_value]
                    if type(old_value) != list:
                        old_value = [old_value]
                    if new_value  == old_value:
                        continue
                    else:
                        merged_value = list(set(new_value + old_value))
                        variant_dict[genome_coor].INFO[key] = copy.deepcopy(merged_value)
    print "number of repeat records: ", num_repeats
#    for value in variant_dict.values():
#        vcf_writer.write_record(value)


    write_to_vcf(args.output, variant_dict)

def write_to_vcf(path_out, v_dict):
    f_out = open(path_out, "w")
    for record in v_dict.values():
        if record.QUAL == None:
            record.QUAL = "."
        if record.FILTER == None:
            record.FILTER = "."

        items = [record.CHROM, str(record.POS), str(record.ID), record.REF, 
                str(record.ALT[0]), record.QUAL, record.FILTER]
        infos = []
        for key in record.INFO:
            this_info = record.INFO[key] 
            if type(this_info) == list:
                if None in this_info:
                    this_info = ['None' if x is None else x for x in this_info]
                infos.append(key + "=" + ",".join(this_info))
            else:
                infos.append(key + "=" + str(this_info))
        items.append(";".join(infos))
        new_line = "\t".join(items)
        f_out.write(new_line + "\n")




if __name__ == "__main__":
    main()

