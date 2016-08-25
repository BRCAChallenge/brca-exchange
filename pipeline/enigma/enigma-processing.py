"""
this scripts takes the enigma raw tsv file and process it to produce the format suitable for 
brcaexchange variant merging pipeline
"""

import glob
import numpy as np
import sys
import re
import argparse
import os
from pprint import pprint as pp
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pygr.seqdb import SequenceFileDB
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import datetime



COLUMNS_TO_SAVE = np.array(["Gene_symbol", #Genomic_Coordinate
                            "Reference_sequence",
                            "HGVS", # change to HGVS_cDNA
                            "BIC_Nomenclature",
                            "Abbrev_AA_change",
                            "URL",
                            "Condition_ID_type",
                            "Condition_ID_value",
                            "Condition_category",
                            "Clinical_significance",
                            "Date_last_evaluated",
                            "Assertion_method",
                            "Assertion_method_citation",
                            "Clinical_significance_citations",
                            "Comment_on_clinical_significance",
                            "Collection_method",
                            "Allele_origin",
                            "ClinVarAccession"]) #"HGVS_protein"

OUTPUT_COLUMNS = [i + "_cDNA" if i == "HGVS" else i for i in COLUMNS_TO_SAVE] + ["HGVS_protein"]
OUTPUT_COLUMNS.insert(1, "Genomic_Coordinate")

HDP = hgvs.dataproviders.uta.connect()
EVM = hgvs.variantmapper.EasyVariantMapper(HDP,
    primary_assembly='GRCh37', alt_aln_method='splign')
HP = hgvs.parser.Parser()


def main():
    raw_files = sorted(glob.glob("raw_files/*batch*tsv"))

    # Luigi moves all output files to pipeline/brca/pipeline-data/pipeline_input.
    # This code is here in case the output directory is needed instead of using the luigi output.
    today = datetime.date.today().isoformat()
    if not os.path.exists("output"):
        os.makedirs("output")
    default_output_file = "output/ENIGMA_last_updated_%s.tsv" % (today)
    # writable_default_output_file = open(default_output_file, "w")

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--writable_output', type=argparse.FileType('w'),
                        help='Opened writable output file for conversion.',
                        default=default_output_file)
    parser.add_argument('-g', '--genome_path',
                        help='Link to hg38.fa.')
    args = parser.parse_args()

    f_out = args.writable_output
    GENOME = SequenceFileDB(args.genome_path)

    f_out.write("\t".join(OUTPUT_COLUMNS) + "\n")
    for filename in raw_files:
        f_in = open(filename, "r")
        # print filename
        for index, line in enumerate(f_in):
            # print index
            if index in [0, 1, 3, 4]:
                continue
            items = np.array(line.rstrip().split("\t"))
            if index == 2:
                columns = np.array([i.replace(" ", "_") for i in items])
                index_to_save = [np.where(columns==i)[0][0] for i in COLUMNS_TO_SAVE]
                column_idx = dict(zip(COLUMNS_TO_SAVE, index_to_save))
                continue
            if len(items) != len(columns):
                continue
            OMIM_id_index = column_idx["Condition_ID_value"]
            items[OMIM_id_index] = convert_OMIM_id(items[OMIM_id_index])
            HGVS_cDNA = (items[column_idx["Reference_sequence"]] +
                         ":" + items[column_idx["HGVS"]])
            try:
                genome_coor, HGVS_p = convert_HGVS(HGVS_cDNA, GENOME)
            except:
                # TODO: better handling of misnamed HGVS string
                genome_coor, HGVS_p = create_None_filler() 
            final_items = list(items[index_to_save])
            final_items.insert(1, genome_coor)
            final_items.append(HGVS_p)
            new_line = "\t".join(list(final_items)) + "\n"
            f_out.write(new_line)
        f_in.close()
    f_out.close()


def convert_OMIM_id(OMIM_id):
    if OMIM_id == "604370":
        return "BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 1; BROVCA1 (" + OMIM_id + ")"
    elif OMIM_id == "612555":
        return "BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (" + OMIM_id + ")"
    else:
        raise Exception("OMIM id not found (" + OMIM_id + ")")


def convert_HGVS(hgvs_c, GENOME):
    chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
        hgvs_c, genome, get_transcript=get_transcript)
    genome_coor = chrom + ":" + str(offset) + ":" + ref + ">" + alt
    HGVS_p = HGVS_cDNA_to_protein(hgvs_c)
    return genome_coor, HGVS_p


def create_None_filler():
    genome_coor = "None:None:None>None"
    HGVS_p = "None"
    return genome_coor, HGVS_p


def HGVS_cDNA_to_protein(hgvs_c):
    hgvs_c = HP.parse_hgvs_variant(hgvs_c)
    hgvs_p = str(EVM.c_to_p(hgvs_c)).split(":")[1]
    return hgvs_p


def get_transcript(name):
    # TODO: Use environment variable
    REFGENE = "../resources/refseq/hg38.BRCA.refGene.txt"
    with open(REFGENE) as infile:
        TRANSCRIPTS = pyhgvs_utils.read_transcripts(infile)
    return TRANSCRIPTS.get(name)


if __name__ == "__main__":
    main()
