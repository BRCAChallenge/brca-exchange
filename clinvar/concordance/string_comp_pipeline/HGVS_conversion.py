import pandas as pd
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pygr.seqdb import SequenceFileDB
import sys
import os
import json
import string_comp
ORIGINAL_FILE = "../BRCA_selectedLabs_only/ClinVarBRCA.selectedLabsOnly.txt"


ERROR = "../BRCA_selectedLabs_only/BRCA.wrong_genome_Coor"
FILE = "../BRCA_selectedLabs_only/BRCA.pre-processed"

GENOME = SequenceFileDB("../reference_files/hg19.fa")
with open('../reference_files/genes.refGene.BRCA.txt') as infile:
    transcripts = pyhgvs_utils.read_transcripts(infile)
def get_transcript(name):
    return transcripts.get(name)

def main():
    check_HGVS_conversion_error()
    #HGVS_conversion()



def HGVS_to_genome_coor(HGVS):
    try:
        chrm, pos, ref, alt = pyhgvs.parse_hgvs_name(
            HGVS, GENOME, get_transcript=get_transcript)
        chrm = chrm[3:]
        pos = str(pos)
        ref = ref.replace("-", "")
        alt = alt.replace("-", "")
        genome_coor = "_".join([chrm, pos, ref, alt])
    except(pyhgvs.InvalidHGVSName, NotImplementedError, ValueError, AssertionError):
        genome_coor = "not translated"
    return genome_coor

def check_HGVS_conversion_error():
    df = pd.read_csv(ORIGINAL_FILE, sep="\t")
    n_diff = 0
    n_real_diff = 0
    for row in df.iterrows():
        original_genome_coor = [row[1].Chrom, row[1].Pos,
                                row[1].Ref.replace("-", ""),
                                row[1].Alt.replace("-", "")]
        try:
            chrm, pos, ref, alt = pyhgvs.parse_hgvs_name(
                row[1].HGVS, GENOME, get_transcript=get_transcript)
            chrm = chrm[3:]
            pos = str(pos)
            converted_genome_coor = [chrm, pos, ref, alt]
        except (pyhgvs.InvalidHGVSName, NotImplementedError, ValueError, AssertionError):
            converted_genome_coor = ["NA"]

        if converted_genome_coor == original_genome_coor:
            continue
        elif "NA" in converted_genome_coor or "None" in original_genome_coor:
            continue
        else:
            n_diff += 1
            if not string_comp.variant_equal(original_genome_coor, converted_genome_coor):
                n_real_diff += 1


    print "total number of rows in file: ", len(df)
    print "number of mismatching result by direct comparison: ", n_diff
    print "number of mismatching result by string comparison: ", n_real_diff



if __name__ == "__main__":
    main()