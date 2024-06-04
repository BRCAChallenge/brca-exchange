#!/usr/bin/env python

import argparse
import csv
import re
import sys

BIOINFO_CODE_ID = "Provisional_Evidence_Code_Bioinfo"
BIOINFO_CODE_DESCR = "Provisional_Evidence_Description_Bioinfo"


NO_CODE = "NO_CODE"
PP3 = "PP3"
BP4_BP7 = "BP4,BP7"
BP4 = "BP4"
BP1_STRONG = "BP1_STRONG"
PVS1_CODE = "PVS1"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="build_final.tsv",
                        help="built_final")
    parser.add_argument("-o", "--output", default="built_with_bioinfo.tsv",
                        help="version of input file with new columns added")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Print debugging info")
    args = parser.parse_args()
    return(args)


def initialize_output_file(input_file, output_filename):
    """
    Create an empty output file with the new columns                        
    """
    new_columns = [BIOINFO_CODE_ID, BIOINFO_CODE_DESCR]
    input_header_row = input_file.fieldnames
    if "change_type" in input_header_row:
        idx = input_header_row.index("change_type")
        output_header_row = input_header_row[:idx] + new_columns \
            + input_header_row[idx:]
    else:
        output_header_row = input_header_row + new_columns
    output_file = csv.DictWriter(open(output_filename,"w"),
                                 fieldnames=output_header_row,
                                 delimiter = '\t')
    output_file.writeheader()
    return(output_file)


def extract_protein_coordinate(variant):
    coordinate = None
    hit = re.search("[0-9]+", variant["Protein_Change"])
    if hit:
        token = variant["Protein_Change"][hit.start():hit.end()]
        pos = int(token)
        print("from", variant["Protein_Change"], "derived", pos)
        return(pos)

def inside_functional_domain(variant):
    inside_domain = False
    pos = extract_protein_coordinate(variant)
    if pos:
        if variant["Gene_Symbol"] == "BRCA1":
            if pos >= 2 and pos <= 99:
                inside_domain = True
            elif pos >= 503 and pos <= 508:
                inside_domain = True
            elif pos >= 607 and pos <= 614:
                inside_domain = True
            elif pos >= 651 and pos <= 656:
                inside_domain = True
            elif pos >= 1391 and pos <= 1424:
                inside_domain = True
            elif pos >= 1650 and pos <= 1863:
                inside_domain = True
        elif variant["Gene_Symbol"] == "BRCA2":
            if pos >= 10 and pos <= 40:
                inside_domain = True
            elif pos >= 1002 and pos <= 1036:
                inside_domain = True
            elif pos >= 1212 and pos <= 1246:
                inside_domain = True
            elif pos >= 1422 and pos <= 1453:
                inside_domain = True
            elif pos >= 1518 and pos <= 1549:
                inside_domain = True
            elif pos >= 1665 and pos <= 1696:
                inside_domain = True
            elif pos >= 1837 and pos <= 1871:
                inside_domain = True
            elif pos >= 1971 and pos <= 2005:
                inside_domain = True
            elif pos >= 2051 and pos <= 2085:
                inside_domain = True
            elif pos >= 2481 and pos <= 3186:
                inside_domain = True
            elif pos >= 3263 and pos <= 3269:
                inside_domain = True
            elif pos >= 3265 and pos <= 3330:
                inside_domain = True
            elif pos >= 3381 and pos <= 3385:
                inside_domain = True
    return(inside_domain)


        
def estimate_bioinfo_code(variant):
    effect = "unknown"
    bioinfo_code = NO_CODE
    if re.search("=\)$", variant["pyhgvs_Protein"]):
        effect = "synonymous_variant"
    elif re.search("[A-Z]+[0-9]+[A-Z]+", variant["Protein_Change"]):
        effect = "missense_variant"
    elif re.search("c\.[0-9]+[+]", variant["pyhgvs_cDNA"]):
        effect = "intron_variant"
    elif re.search("c\.[0-9]+[-]", variant["pyhgvs_cDNA"]):
        effect = "intron_variant"
    print("variant", variant["pyhgvs_cDNA"], "protein change", variant["Protein_Change"], variant["pyhgvs_Protein"], "effect", effect)
    if  variant["result_spliceai"] == "-":
        splicing_effect = False
        no_splicing_effect = True
    else:
        splicing_effect = (float(variant["result_spliceai"]) > 0.2)
        no_splicing_effect = (float(variant["result_spliceai"]) < 0.1)
    if variant["Gene_Symbol"] == "BRCA1":
        if variant["BayesDel_nsfp33a_noAF"] == "-":
            protein_effect = False
            no_protein_effect = True
        elif float(variant["BayesDel_nsfp33a_noAF"]) > 0.28:
            protein_effect = True
            no_prptein_effect = False
        elif float(variant["BayesDel_nsfp33a_noAF"]) < 0.15:
            protein_effect = False
            no_protein_effect = True
        else:
            protein_effect = False
            no_protein_effect = False
    if variant["Gene_Symbol"] == "BRCA2":
        if variant["BayesDel_nsfp33a_noAF"] == "-":
            protein_effect = False
            no_protein_effect = True
        elif float(variant["BayesDel_nsfp33a_noAF"]) > 0.30:
            protein_effect = True
            no_prptein_effect = False
        elif float(variant["BayesDel_nsfp33a_noAF"]) < 0.18:
            protein_effect = False
            no_protein_effect = True
        else:
            protein_effect = False
            no_protein_effect = False
    inside_domain = inside_functional_domain(variant)
    print("effect", effect, "splicing effect", splicing_effect, "inside domain", inside_domain)
    if effect == "synonymous_variant":
        if splicing_effect:
            bioinfo_code = PP3
        elif inside_domain:
            bioinfo_code = BP4_BP7
        else:
            bioinfo_code = BP1_STRONG
    elif effect == "intron_variant":
        if splicing_effect:
            bioinfo_code = PP3
        else:
            bioinfo_code = BP4
    elif effect == "missense_variant":
        if splicing_effect:
            bioinfo_code = PP3
        elif no_splicing_effect:
            if not inside_domain:
                bioinfo_code = BP1_STRONG
            elif protein_effect:
                bioinfo_code = PP3
            elif no_protein_effect:
                bioinfo_code = BP4
        else:
            if inside_domain and protein_effect:
                bioinfo_code = PP3
    return(bioinfo_code)


def apply_pvs1_code(variant):
    pvs1_code = NO_CODE
    protein_hgvs = variant["pyhgvs_Protein"]
    stop_added = re.search("Ter", protein_hgvs)
    if stop_added:
        pvs1_code = PVS1_CODE
    return(pvs1_code)

    
def main():
    csv.field_size_limit(sys.maxsize)
    args = parse_args()
    with open(args.input, 'r') as input_fp:
        input_reader = csv.DictReader(input_fp, delimiter = "\t")
        writer = initialize_output_file(input_reader, args.output)
        for variant in input_reader:
            #variant[BIOINFO_CODE_ID] = estimate_bioinfo_code(variant, debug=args.debug)
            #pvs1_code = apply_pvs1_code(variant)
            variant[BIOINFO_CODE_ID] = ""
            variant[BIOINFO_CODE_DESCR] = ""
            writer.writerow(variant)

if __name__ == "__main__":
    main()
