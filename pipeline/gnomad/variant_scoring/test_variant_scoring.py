#!/usr/bin/env python
# coding: utf-8

import itertools
import argparse
import csv
from collections import OrderedDict
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="built_with_popfreq.tsv",
                        help="Input file with new estimated popfreq")
    parser.add_argument("-t", "--test", default="test.tsv",
                        help="Input file with old estimted popfreq")
    args = parser.parse_args()
    return(args)

def field_defined(field):
    """
    Return a binary value indicating whether or not this field 
    has a defined value
    """
    return(field != "-")

def key_fields(id, variant, evidence_code, faf95, allele_count, description):
    data_this_variant = {
        "id" : id,
        "hgvs" : variant["pyhgvs_cDNA"],
        "evidence_code" : evidence_code,
        "faf95": faf95,
        "allele_count": allele_count,
        "description": description,
        "final_code" : variant["Provisional_evidence_code_popfreq"]
        }
    return(data_this_variant)


def parse_gnomad_annotations(input_file):
    v2_new = dict()
    v3_new = dict()
    for variant in input_file:
        if field_defined(variant["Variant_id_GnomAD"]):
            v2_id = variant["Variant_id_GnomAD"]
            v2_new[v2_id] = key_fields(v2_id, variant,
                                       variant["Provisional_code_GnomAD"],
                                       variant["faf95_popmax_exome_GnomAD"],
                                       variant["Allele_count_exome_GnomAD"],
                                       variant["Provisional_code_description_GnomAD"])
        if field_defined(variant["Variant_id_GnomADv3"]):
            v3_id = variant["Variant_id_GnomADv3"]
            v3_new[v3_id] = key_fields(v3_id, variant,
                                       variant["Provisional_code_GnomADv3"],
                                       variant["faf95_popmax_genome_GnomADv3"],
                                       variant["Allele_count_genome_GnomADv3"],
                                       variant["Provisional_code_description_GnomADv3"])
    return(v2_new, v3_new)

def print_header():
    print("\t".join(["gnomAD_id", "gnomAD_version", "hgvs",
                     "evidence_code_old", "faf95_old", "allele_freq_old",
                     "evidence_code_new", "faf95_new", "allele_freq_new",
                     "descr_new",
                     "final_code_old", "final_code_new"]))

def compare(old_data, version, var_name, new_data):
    if var_name in new_data:
        print("\t".join([var_name, version, new_data[var_name]["hgvs"],
                         old_data["evidence_code"],
                         old_data["AF_popmax"],
                         new_data[var_name]["evidence_code"],
                         new_data[var_name]["faf95"],
                         new_data[var_name]["allele_count"]
                         new_data[var_name]["description"],
                         old_data["final_code"],
                         new_data[var_name]["final_code"]]))
    
    
def print_comparison(entry, v2, v3):
    version = entry["src"]
    var_name = entry["var_name"]
    if version == "v2":
        compare(entry, version, entry["var_name_hg19"], v2)
    elif version == "v3":
        compare(entry, version, var_name, v3)
        
    
def main():
    args = parse_args()
    input_file = csv.DictReader(open(args.input), delimiter = "\t")
    (v2_new, v3_new) = parse_gnomad_annotations(input_file)
    test_file = csv.DictReader(open(args.test), delimiter = "\t")
    print_header()
    for entry in test_file:
        print_comparison(entry, v2_new, v3_new)
    

if __name__ == "__main__":
    main()
