"""
this script takes allele frequency input from exac, 1000g and esp and add allele frequency to
BRCA table
"""
import pandas as pd
import string_comp
import sys
import pickle
from pprint import pprint as pp

PATH = "/Users/Molly/Desktop/BRCA Research/invitae_paper/"

IN_FILE = PATH + "BRCA_selectedLabs_only/BRCA.transposed_by_uniq_id.add_concordance"
OUT_FILE = PATH + "BRCA_selectedLabs_only/table_BRCA_ucsc_string_comparison.tsv"
DF = pd.read_csv(IN_FILE, sep="\t", dtype=str)

EXAC = PATH + "allele_frequency/exac_af.txt"
ONEKG = PATH + "allele_frequency/1000g_af.txt"
ESP_AA = PATH + "allele_frequency/esp_aa_af.txt"
ESP_EA = PATH + "allele_frequency/esp_ea_af.txt"

SOURCES = {"esp_aa": ESP_AA, "esp_ea": ESP_EA, "exac_af": EXAC, "1000g_af": ONEKG}
SERIALIZED = PATH + "BRCA_selectedLabs_only/serialized_uniq_id_mapping"



def main():
    #check_af_source()
    uniq_ids = list(DF.uniq_id)
    new_columns = list(DF.columns)
    for source_name, source_file in SOURCES.iteritems():
        new_columns.append(source_name)
        DF[source_name] = pd.Series(get_allele_frequency(source_file, uniq_ids))

    new_df = DF[new_columns]
    new_df.to_csv(OUT_FILE, sep="\t", dtype=str, index=False)


def get_allele_frequency(source, uniq_ids):
    df = pd.read_csv(source, sep="\t", header=None, dtype=str, names=["genome", "af"])
    af_dict = dict(zip(df.genome, df.af))
    af_list = []
    for uniq_id in uniq_ids:
        this_af = "x"
        print len(af_list)
        if "NM" in uniq_id:
            af_list.append(this_af)
            continue

        if "," in uniq_id:
            this_af = get_af_of_merged_variant(uniq_id, af_dict)
            af_list.append(this_af)
            continue

        if uniq_id in af_dict:
            this_af = af_dict[uniq_id]
            af_list.append(this_af)
            continue
        else:
            for each_genome_coor in af_dict:
                v1 = uniq_id.split("_")
                v2 = each_genome_coor.split("_")
                if string_comp.variant_equal(v1, v2):
                    this_af = af_dict[each_genome_coor]
                    break
        af_list.append(this_af)
    assert(len(af_list) == len(DF))
    return af_list

def get_af_of_merged_variant(uniq_id, af_dict):
    variants = uniq_id.split(",")
    for v in variants:
        if v in af_dict:
            return af_dict[v]
    for each_genome_coor in af_dict:
        v1 = variants[0].split("_")
        v2 = each_genome_coor.split("_")
        if string_comp.variant_equal(v1, v2):
            return af_dict[each_genome_coor]
    return "x"




def check_af_source():
    """
    check two things:
    1) if each allele frequency source has equivalent variants
    2) if each source file has wrong genomic coordinates
    """
    for source in SOURCES.values():
        print("checking " + source)
        df = pd.read_csv(source, sep=" ", header=None, names = ["genome", "af"])
        # check genomic coordinate correctness
        for v in df.genome:
            v1 = v.split("_")
            if not string_comp.ref_correct(v1):
                raise Exception("wrong genomic coordinate in {0}: {1}".format(v, source))
        # check existence of equivalent variants
        v_list = list(df.genome)
        while len(v_list) > 1:
            if len(v_list)%10 == 0:
                print len(v_list)
            first = v_list.pop(0)
            for each_rest in v_list:
                if string_comp.variant_equal(first.split("_"), each_rest.split("_")):
                    raise Exception("equivalent variant detected in {0}".format(source))

        print "finished checking {0}".format(source)
        print "all good: no wrong genomic coordinates or equivalent variants"
    print "Done"



if __name__ == "__main__":
    main()

