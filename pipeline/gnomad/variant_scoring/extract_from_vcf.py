#!/usr/bin/env python                                                                                                   

"""                                                                                                                      
Description:                                                                                                            
    Takes in a gnomad table and converts it to vcf format.                                                               
"""

import argparse
from pysam import VariantFile

def parse_args():
    parser = argparse.ArgumentParser(description="Parse a gnomAD VCF for the BRCA Exchange pipeline")
    parser.add_argument('-i', '--input', help="Input VCF file, from gnomAD")
    parser.add_argument('-v', '--verbose', action='count', default=False, help='determines logging')
    options = parser.parse_args()
    return options

def main():
    args = parse_args()
    vcf_in = VariantFile(args.input)
    print("ID\tAC_joint\tAF_joint\tAN_joint\tfafmax_faf95_max_joint" \
          + "\tfafmax_faf95_max_gen_anc_joint" \
          + "\tjoint_filters\tgenome_filters\texome_filters" \
          + "\tnot_called_in_genome\tnot_called_in_exome")
    for record in vcf_in.fetch():
        all_filters = ""
        delim = ""
        for this_filter in record.filter.keys():
            all_filters = all_filters + delim + this_filter
            delim = ","
        alt_idx = 0
        for this_alt in record.alts:
            if "fafmax_faf95_max_joint" in record.info:
                fafmax_faf95 = str(record.info["fafmax_faf95_max_joint"][alt_idx])
            else:
                fafmax_faf95 = "NA"
            if "fafmax_faf95_max_gen_anc_joint" in record.info:
                fafmax_max_gen_anc = record.info["fafmax_faf95_max_gen_anc_joint"][alt_idx]
            else:
                fafmax_max_gen_anc = "NA"
            if "genomes_filters" in record.info:
                genomes_filters = record.info["genomes_filters"]
            else:
                genomes_filters = ""
            if "exomes_filters" in record.info:
                exomes_filters = record.info["exomes_filters"]
            else:
                exomes_filters = ""
            print(("%s-%s-%s-%s" + "\t%d" + "\t%f\t%d" + "\t%s\t%s\t%s"
                   + "\t%s\t%s" + "\t%s" + "\t%s")
                  % (record.chrom, record.pos, record.ref, this_alt,
                     record.info["AC_joint"][alt_idx],
                     record.info["AF_joint"][alt_idx],record.info["AN_joint"],
                     fafmax_faf95, fafmax_max_gen_anc, all_filters,
                     genomes_filters, exomes_filters,
                     str(record.info["not_called_in_genomes"]),
                     str(record.info["not_called_in_exomes"])))
            alt_idx += 1
    
if __name__ == "__main__":
    main()
