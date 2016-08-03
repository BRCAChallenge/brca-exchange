"""
this script compares the concordance of variant submitters by comparing variants with three methods:
1. by variant ID
2. by genomic coordinate converted to DNA level (string comparison)
"""

import operator
import sys
import argparse
import string_comp
import transpose_variant as tv
import pandas as pd


CLINSIGS = ["Benign", "Likely benign", "Likely pathogenic", 
            "Pathogenic", "Uncertain significance"]
SUBMITTERS = ["Invitae", "GeneDx", "Ambry", "Emory", "Children", 
              "SCRP", "Counsyl"]
IN_FILE = "../BRCA_selectedLabs_only/BRCA.pre-processed.add_uniq_id"

parser = argparse.ArgumentParser(description='choose verbosity level')
parser.add_argument("-v", "--verbose",
                    help="increase output verbosity", action="store_true")
args = parser.parse_args()

def main():
    analysis_by_string_genomic_coor()

def analysis_by_string_genomic_coor():
    all_dict = {}

    print "number of equivalent variant within one submitter (discordant submissions)"
    for submitter in SUBMITTERS:
        this_dict = separate_submitter(submitter, True)
        all_dict[submitter] = this_dict

    print "\nnumber of unique variants from each submitter:"
    for key, value in all_dict.iteritems():
        print key, len(value)

    print "\nsubmitter1 and submitter2: common variants (discordant variants)"
    calculate_concordance(all_dict)

def calculate_concordance(all_dict):
    compared_pairs = []
    for submitter in all_dict.keys():
        rest_submitter_list = all_dict.keys()
        rest_submitter_list.remove(submitter)
        for rest_submitter in rest_submitter_list:
            if sorted([submitter, rest_submitter]) in compared_pairs:
                continue
            compare_two_dict(submitter, rest_submitter, all_dict)
            compared_pairs.append(sorted([submitter, rest_submitter]))
    
def compare_two_dict(submitter1, submitter2, all_dict):
    df = pd.read_csv(IN_FILE, sep="\t", dtype=str)
    d1 = all_dict[submitter1]
    d2 = all_dict[submitter2]
    common = set(d1.keys()) & set(d2.keys())
    disagreement = set()
    discordant = set()
    for key in common:
        if d1[key] != d2[key]:
            disagreement.add(key)
        if is_discordant(d1[key], d2[key]):
            discordant.add(key)
    print "{0} and {1}: {2} ({3})".format(submitter1, submitter2,
                                          len(common), len(discordant))
    if len(discordant) != 0 and args.verbose:
        print "-------------details----------------------"
        for each_id in discordant:
            rows_with_id = df[df.uniq_id == each_id]
            print "id: ", each_id
            for index, row in rows_with_id.iterrows():
                print ",".join([row.HGVS, row.Submitter, row.ClinicalSignificance, str(row.VariantID)])
        print "------------------------------------------"

def is_discordant(term1, term2):
    """assess concordinance based on invitae poster criteria"""
    # if (term1 not in CLINSIGS) or (term2 not in CLINSIGS):
    #     raise Exception("bad clincical description input")
   
    if "PATHO" in (term1.upper() + term2.upper()):
        if ("PATHO" in term1.upper()) and ("PATHO" in term2.upper()):
            return False
        else:
            return True
    else:
        return False

def separate_submitter(submitter_abbrev, string_level):
    """save variant id from one submitter in one dictionary
    !!! if one submitter has multiple submission of the same variant,
    only the last submision is recorded
    """
    df = pd.read_csv(IN_FILE, sep="\t", dtype=str)
    submitter_dict = {}
    repeat_id = {}
    n_repeat_ids = 0
    for index, row in df.iterrows():
        if string_level:
            this_id = row.uniq_id
        else:
            this_id = row.VariantID
        if submitter_abbrev in row.Submitter:
            if this_id not in repeat_id:
                repeat_id[this_id] = [row.ClinicalSignificance]
            else:
                repeat_id[this_id].append(row.ClinicalSignificance)
                n_repeat_ids += 1
            submitter_dict[this_id] = row.ClinicalSignificance

    n_discordance = 0
    for each_id, pathos in repeat_id.iteritems():
        if len(pathos) > 1:
            if not tv.decide_concordance(set(pathos)):
                n_discordance += 1

    print "{0}: {1} ({2})".format(submitter_abbrev, n_repeat_ids, n_discordance)

    return submitter_dict

if __name__=="__main__":
    main()

