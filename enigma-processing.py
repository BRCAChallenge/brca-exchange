"""
this scripts takes the enigma variant list and do various things with it, inlcuding:
    1. convert HGVS to genomic coordinate
"""

import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import pyhgvs



def main():
    enigma_file = open("data/Enigma_variants/enigma-database.tsv", "r")
    line_num = 0
    for line in enigma_file:
        line_num += 1
        if line_num == 1:
            pass 
        else:
            items = line.strip().split("\t")
            transcript_id = items[1]
            HGVS = items[2]
            hgvs_c = transcript_id + ":" + HGVS
            genome = convert_hgvs_to_genome(hgvs_c)
            items.append(genome)

def convert_hgvs_to_genome(hgvs_c):
    hp = hgvs.parser.Parser()
    print hgvs_c
    hgvs_c = hp.parse_hgvs_variant(hgvs_c)
    hdp = hgvs.dataproviders.uta.connect()
    evm = hgvs.variantmapper.EasyVariantMapper(hdp,
        primary_assembly='GRCh37', alt_aln_method='splign')
    hgvs_g = evm.c_to_g(hgvs_c)
    hgvs_p = evm.c_to_p(hgvs_c)
    hgvs_n = evm.c_to_n(hgvs_c)
    print hgvs_n
    print " "





if __name__ == "__main__":
    main()
