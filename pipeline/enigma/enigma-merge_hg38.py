'''
This script takes processed ENIGMA files and concatenates them.
A list of files can be taken as a command line argument or the
default using the output file naming convention from the processing scripts
'''

import glob
import numpy as np
import sys
import re
import argparse
import os
import datetime

COLUMNS = np.array(["Gene_symbol",
                    "Genomic_Coordinate",
                    "Reference_sequence",
                    "HGVS_cDNA",
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
                    "ClinVarAccession",
                    "HGVS_protein"])


def main():
    default_input_files = sorted(glob.glob("output/*last_updated*hg38*tsv"))
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_folder', help='Link to folder with all processed Enigma files to merge')
    parser.add_argument('-o', '--output_folder', help='Folder to contain output file')
    parser.add_argument('-f', '--input_files', nargs="*",
                        help='Links to processed ENIGMA files to be merged if input folder is not provided.',
                        default=default_input_files)
    args = parser.parse_args()
    if args.input_folder is not None:
        input_files = sorted(glob.glob(args.input_folder + "*last_updated*hg38*tsv"))
    else:
        input_files = (args.input_files)
    output = open(args.output_folder + "ENIGMA_combined_hg38.tsv", "w")
    output.write("\t".join(COLUMNS) + "\n")

    print "input_files", input_files
    for filename in input_files:
        print "combining", filename
        f_in = open(filename, "r")
        for index, line in enumerate(f_in):
            items = np.array(line.rstrip().split("\t"))
            if index == 0:
                index_to_save = [np.where(COLUMNS == i)[0][0] for i in COLUMNS]
                column_idx = dict(zip(COLUMNS, index_to_save))
                continue
            if len(items) != len(COLUMNS):
                continue
            final_items = list(items[index_to_save])
            new_line = "\t".join(list(final_items)) + "\n"
            output.write(new_line)
        f_in.close()
    output.close()

if __name__ == "__main__":
    main()
