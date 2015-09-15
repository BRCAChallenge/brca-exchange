"""
this scripts takes the enigma variant list and do following things:
    1.  append HGVS_cDNA and Genomic Coordinate column
    2. change "comment on classification" column to parsed succinct words
"""

import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
import re
import sys
from pygr.seqdb import SequenceFileDB

def main():
    processed_file = open("data/enigma_database_v3.tsv", "w")
    line_num = 0
    for line in sys.stdin:
        line_num += 1
        print line_num
        if line_num == 1:
            items = line.strip().split("\t")
            items[2] = "HGVS_cDNA"
            items.append("HGVS_protein")
            items.append("Genomic Coordinate") 
            items[13] = "Classification method"
        else:
            items = line.strip().split("\t")
            transcript_id = items[1]
            HGVS = items[2]
            hgvs_c = transcript_id + ":" + HGVS
            converted_hgvs = convert_hgvs(hgvs_c)
            items += converted_hgvs
            items[13] = parse_comment(items[13])
        new_line = "\t".join(items) + "\n"
        processed_file.write(new_line)

def parse_comment(s):
    text = ""
    if "IARC" in s:
        text = re.search("posterior probability = (([0-9.]+)|1)", s)
        text = text.group(0).strip(".")
    if "IARC" in s and "frequency" in s:
        text += " and frequency >1%"
    elif "frequency" in s and "IARC" not in s:
        text = "frequency > 1%"
    return text


def convert_hgvs(hgvs_c):
    hp = hgvs.parser.Parser()
    hgvs_c = hp.parse_hgvs_variant(hgvs_c)
    hdp = hgvs.dataproviders.uta.connect()
    evm = hgvs.variantmapper.EasyVariantMapper(hdp,
        primary_assembly='GRCh37', alt_aln_method='splign')
    #hgvs_g = evm.c_to_g(hgvs_c)
    hgvs_p = evm.c_to_p(hgvs_c)
    hgvs_p_no_transcript = str(hgvs_p).split(":")[1]
    genome_coordinate = get_genome_coor(str(hgvs_c))
    return [hgvs_p_no_transcript, genome_coordinate]

def get_genome_coor(hgvs_c):
    genome = SequenceFileDB('data/hg19.fa')
    refGene = "/Users/Molly/Desktop/web-dev/hgvs_counsyl/hgvs/pyhgvs/data/genes.refGene"
    with open(refGene) as infile:
        transcripts = pyhgvs_utils.read_transcripts(infile)

    def get_transcript(name):
        return transcripts.get(name)
    
    chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
        hgvs_c, genome, get_transcript=get_transcript)
    return chrom + ":" + str(offset) + ":" + ref + ">" + alt



def parse_hgvs_g(hgvs_g):
    text = re.search("NC_0000(1[37])\.10:g\.([0-9]+)([AGTC]>[AGTC])", hgvs_g)
    chrom = text.group(1)
    position = text.group(2)
    change = text.group(3)
    return "chr" + chrom + ":" + position + ":" + change

if __name__ == "__main__":
    main()
