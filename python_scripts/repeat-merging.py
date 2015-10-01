"""
this script takes a vcf file and collapes repeatative variant rows
"""

import vcf
import argparse
import copy

def main():
    args = arg_parse()
    vcf_reader = vcf.Reader(open(args.input, "r"), strict_whitespace=True)
    vcf_writer = vcf.Writer(open(args.output, 'w'), vcf_reader)
    variant_dict = {}
    num_repeats = 0
    for record in vcf_reader:
        genome_coor = "chr{0}:{1}:{2}>{3}".format(
            record.CHROM, str(record.POS), record.REF, record.ALT[0])
        if genome_coor not in variant_dict.keys():
            variant_dict[genome_coor] = copy.deepcopy(record)
        else:
            num_repeats += 1
#           print "before -----------"
#           print variant_dict[genome_coor]
#           print variant_dict[genome_coor].INFO
#           print record
#           print record.INFO
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
#           print "after -----------"
#           print variant_dict[genome_coor]
#           print variant_dict[genome_coor].INFO
#           print record
#           print record.INFO
    print "number of repeat records: ", num_repeats
    for value in variant_dict.values():
        vcf_writer.write_record(value)



def arg_parse():
    parser = argparse.ArgumentParser() 
    parser.add_argument("-i", "--input") 
    parser.add_argument("-o", "--output") 
    args = parser.parse_args() 
    return args 



if __name__ == "__main__":
    main()

