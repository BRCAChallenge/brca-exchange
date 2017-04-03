"""
NOTES:

The input file must have the appropriate column and formatting adjustments as outlined in ./README.md
It will not run on ENIGMA data prior to 10-18-2016 without making the aforementioned changes.

This script takes a raw enigma tsv file with clinvar accessions and processes it to produce a file that's suitable for
the brcaexchange pipeline.
"""

import glob
import numpy as np
import sys
import re
import argparse
import os
from Bio.SeqUtils import seq1
from pprint import pprint as pp
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pygr.seqdb import SequenceFileDB
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import datetime


COLUMNS_TO_SAVE = np.array(["Gene_symbol",
                            "Reference_sequence",
                            "HGVS",
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
                            "ClinVarAccession"])

OUTPUT_COLUMNS = [i + "_cDNA" if i == "HGVS" else i for i in COLUMNS_TO_SAVE] + ["HGVS_protein"]
OUTPUT_COLUMNS.insert(1, "Genomic_Coordinate")
HDP = hgvs.dataproviders.uta.connect()
EVM = hgvs.variantmapper.EasyVariantMapper(HDP, primary_assembly='GRCh37', alt_aln_method='splign')
HP = hgvs.parser.Parser()
REFGENE = None


def main():
    global REFGENE

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--readable_input', 
                        help='readable input file for conversion.')
    parser.add_argument('-o', '--writable_output', 
                        help='writable output file for conversion.')
    parser.add_argument('-g', '--genome_path', help='Link to hg38.fa.')
    parser.add_argument('-r', '--reference_genome', default='./hg38.BRCA.refGene.txt',
                        help='Link to hg38.BRCA.refgene.txt.')

    args = parser.parse_args()
    GENOME = SequenceFileDB(args.genome_path)
    REFGENE = args.reference_genome

    f_in = open(args.readable_input, "r")
    f_out = open(args.writable_output, "w")
    f_out.write("\t".join(OUTPUT_COLUMNS) + "\n")
    for index, line in enumerate(f_in):
        # 
        # Clean the line by removing leading or trailing spaces adjacent to tabs.  
        #
        line = re.sub("( )*\t( )*", "\t", line)
        items = np.array(line.rstrip().split("\t"))
        if index == 0:
            # Handle column names
            columns = np.array([i.replace(" ", "_") for i in items])
            index_to_save = [np.where(columns == i)[0][0] for i in COLUMNS_TO_SAVE]
            column_idx = dict(zip(COLUMNS_TO_SAVE, index_to_save))
            continue
        #
        # In the date last evaluated field, delete the time last evaluated if provided.
        #
        date_last_evaluated_idx = column_idx["Date_last_evaluated"]
        items[date_last_evaluated_idx] = items[date_last_evaluated_idx].split(' ')[0]
        OMIM_id_index = column_idx["Condition_ID_value"]
        items[OMIM_id_index] = convert_OMIM_id(items[OMIM_id_index])
        items[column_idx["HGVS"]] = cleanup_HGVS(items[column_idx["Reference_sequence"]],
                                 items[column_idx["HGVS"]], HP, EVM)
        HGVS_cDNA = items[column_idx["Reference_sequence"]] + ":" + items[column_idx["HGVS"]]
        print items[column_idx["Reference_sequence"]], items[column_idx["HGVS"]], HGVS_cDNA
        try:
            genome_coor, HGVS_p = convert_HGVS(HGVS_cDNA, GENOME)
        except:
            if (items[column_idx["HGVS"]]).find(";") > -1:
                genome_coor, HGVS_p = create_None_filler()
        aa_abrev_index = column_idx["Abbrev_AA_change"]
        if HGVS_p not in ["p.?", "p.(=)", "None"]:
            if items[aa_abrev_index] == '':
                items[aa_abrev_index] = HGVS_p_to_AA_abrev(HGVS_p)
        final_items = list(items[index_to_save])
        final_items.insert(1, genome_coor)
        final_items.append(HGVS_p)
        new_line = "\t".join(list(final_items)) + "\n"
        f_out.write(new_line)
    f_in.close()
    f_out.close()

def cleanup_HGVS(reference_sequence, cdna_hgvs, hp, evm):
    """The pyhgvs library breaks on HGVS strings that do not specify a ref or alt allele.  The hgvs
    library is robust to these HGVS strings, but does not map to hg38.  Here. use the HGVS library
    to convert the HGVS cDNA string to a HGVS genomic string (against GRCh37), and convert that genomic
    string back to a HGVS cDNA string.  The resulting string can be parsed by the pyhgvs library more
    reliably.  But, it has a strange error with certain insertions.  So don't do this with insertions..."""
    if re.search("ins", cdna_hgvs) and not re.search("del", cdna_hgvs):
        return cdna_hgvs
    else: 
        initial_cdna_hgvs = hp.parse_hgvs_variant(reference_sequence + ":" + cdna_hgvs)
        genomic_hgvs = evm.c_to_g(initial_cdna_hgvs)
        refined_cdna_hgvs = evm.g_to_c(genomic_hgvs, reference_sequence)
        printable_cdna_hgvs = "{xx}".format(xx=refined_cdna_hgvs)
        cleaned_nucleotide_cdna = re.sub(reference_sequence + ":", "", printable_cdna_hgvs)
        return cleaned_nucleotide_cdna



def convert_OMIM_id(OMIM_id):
    if OMIM_id == "604370":
        return "BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 1; BROVCA1 (" + OMIM_id + ")"
    elif OMIM_id == "612555":
        return "BREAST-OVARIAN CANCER, FAMILIAL, SUSCEPTIBILITY TO, 2; BROVCA2 (" + OMIM_id + ")"
    else:
        raise Exception("OMIM id not found (" + OMIM_id + ")")


def convert_HGVS(hgvs_c, GENOME):
    chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
        hgvs_c, GENOME, get_transcript=get_transcript)
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


def HGVS_p_to_AA_abrev(HGVS_p):
    mut = HGVS_p.split("p.(")[1].split(")")[0]
    num = re.findall(r"[+-]?\d+(?:\.\d+)?", mut)
    first_part = mut[0:(len(num[0])+6)]
    frameshift = re.search('fs\w+', mut)
    if frameshift is not None:
        fs = frameshift.group(0)
        aa3 = seq1(mut[-(len(num[1])+3):-len(num[1])])
        fs_part = fs[0:2]+aa3+num[1]
    else:
        fs_part = ''
    aa1 = seq1(first_part[0:3])
    aa2 = seq1(first_part[-3:])
    aa_abbrev = aa1+num[0]+aa2+fs_part
    return aa_abbrev


def get_transcript(name):
    global REFGENE
    if REFGENE is None:
        sys.exit("No reference genome was provided. Try to locate hg38.BRCA.refGene.txt.")
    with open(REFGENE) as infile:
        TRANSCRIPTS = pyhgvs_utils.read_transcripts(infile)
    return TRANSCRIPTS.get(name)


if __name__ == "__main__":
    main()
