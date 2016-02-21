"""
This script expands the numerous subfields of VEP consequence in exac vcf file. 
each of the CSQ subfileds become a new INFO tag starting with CSQ, and the whole VCF file is rearranged
"""

import argparse


VEP = "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|HGNC_ID|BIOTYPE|CANONICAL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|EXON|INTRON|DOMAINS|HGVSc|HGVSp|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EUR_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF"

VEP_fields = VEP.split("|")

def arg_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    return args

def main():
    args = arg_parse()
    f_in = open(args.input, "r")
    f_out = open(args.output, "w")
    #f_out = open("../data/allVcf/exac.VEP_expanded.vcf","w")
    line_num = 0
    for line in f_in:
        line_num += 1
        if line_num <= 193 or line_num == 195:
            f_out.write(line)
        elif line_num == 194:
            for VEP_field in VEP_fields:
                new_line = "##INFO=<ID=CSQ_{0},Number=.,Type=String,Description=\"Consequence type as predicted by VEP.\">\n".format(VEP_field)
                f_out.write(new_line)
        else:
            line_items = line.strip().split("\t")
            INFO = line_items[-1].split(";")
            VEP_string = INFO[-1]
            INFO.pop()

            VEP_list = VEP_string[4:].split(",")
            VEP_list_list = []
            for VEP in VEP_list:
                VEP_results = VEP.split("|")
                VEP_list_list.append(VEP_results)
            merged_list = VEP_list_list[0]
            for alist in VEP_list_list:
                for i in range(49):
                    if alist[i] not in merged_list[i]:
                        merged_list[i] = merged_list[i] + "," + alist[i]
            for i in range(49):
                INFO.append("CSQ_" + VEP_fields[i] + "=" + merged_list[i])

            INFO_string = ";".join(INFO)
            line_items[-1] = INFO_string
            new_line = "\t".join(line_items) + "\n"
            f_out.write(new_line)

if __name__ == "__main__":
    main()

