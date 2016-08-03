"""
this file add a uniq id to each variant row so that equivalent variants have same id
while non-equivalent variants have different id
"""
import pandas as pd
import string_comp

EV = "../BRCA_selectedLabs_only/equivalent_variants.txt"
IN_FILE = "../BRCA_selectedLabs_only/BRCA.pre-processed"

def main():
    #find_equivalent_variant()
    #uniq_ids = check_equivalent_variant_uniqness()
    #add_uniq_id_to_data(uniq_ids)
    get_genome_coor_to_HGVS_dict()

def find_equivalent_variant():
    df = pd.read_csv(IN_FILE, sep="\t", dtype="str")
    set_of_genome_coor = set()
    for index, row in df.iterrows():
        if row.genome_coor != "not translated":
            set_of_genome_coor.add(row.genome_coor)

    uniq_variants = {}
    print "set of of genomic coordinates: ",len(set_of_genome_coor)
    for v in set_of_genome_coor:
        variant_exist = False
        for existing_v in uniq_variants:
            if v == existing_v:
                continue
            else:
                v1 = v.split("_")
                v2 = existing_v.split("_")
                if string_comp.variant_equal(v1, v2):
                    variant_exist = True
                    uniq_variants[existing_v].add(v)
                    print "these two variants are equivlaent", v1, v2
        if not variant_exist:
            uniq_variants[v] = set([v])
    f = open(EV, "w")
    for key, value in uniq_variants.iteritems():
        if len(value) > 1:
            f.write(",".join(list(value)) + "\n")
    f.close()

def add_uniq_id_to_data(ev):
    df = pd.read_csv(IN_FILE, sep="\t", dtype=str)
    uniq_ids = []
    for index, row in df.iterrows():
        if row.genome_coor == "not translated":
            uniq_ids.append(row.HGVS)
        else:
            uniq_ids.append(get_uniq_id(row.genome_coor, ev))
    df['uniq_id'] = pd.Series(uniq_ids)
    df.to_csv(IN_FILE + ".add_uniq_id", sep="\t", dtype=str, index=False)


def get_uniq_id(genome_coor, ev):
    uniq_id = genome_coor
    for each_id in ev:
        if genome_coor in each_id.split(","):
            uniq_id = each_id
    return uniq_id

def get_genome_coor_to_HGVS_dict():
    df = pd.read_csv(IN_FILE, sep="\t", dtype="str")
    x = dict(zip(df.genome_coor, df.HGVS))
    print x['not translated']


def check_equivalent_variant_uniqness():
    """make sure all equivalent variants only occured once
       in the 'equivalent_variants' file"""
    f = open(EV, "r")
    variants = list()
    uniq_id = list()
    for line in f:
        variants += line.strip().split(",")
        uniq_id.append(line.strip())
    f.close()
    if len(variants) != len(set(variants)):
        raise Exception("equivalent_variant file not good")
    print "number of equivalent variants: ", len(variants)

    return uniq_id

if __name__ == "__main__":
    main()
