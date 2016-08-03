"""
improved algorithm for comparing two variants: at DNA string level
TODO: this scripts runs really slow, needs performance improvement
"""

import glob

# remember the trailing "/" otherwise path won't work
PATH = "/hive/groups//brca/phase1/data/cutoff_vcf/"

chr13 = open("brca2.txt", "r")
BRCA2 = chr13.read()
chr17 = open("brca1.txt", "r")
BRCA1 = chr17.read()
BRCA2_START = 32800000
BRCA1_START = 41100000


def main():
    databases = get_databases(PATH)
    for key in databases.keys():
        print_variant_size(databases[key], key)
        for key2 in databases.keys():
            if key != key2:
                get_overlap(key, key2, databases)
        print ""

def get_databases(path):
    db = {}
    files = glob.glob(path + "*.vcf")
    for file in files:
        db_name = file.split(".")[0].split("/")[-1].split(".")[0]
        db[db_name] = get_unique_variants(file)
    return db

def print_variant_size(database, database_name):
    print "unique variants in %s:%d" %(database_name, len(database))

def get_overlap(name_db1, name_db2, database):
    num_overlap = 0
    num_variant = 0
    for variant in database[name_db1]:
        num_variant += 1
        if check_variant_exist(variant, database[name_db2]) == True:
            num_overlap += 1
    print "overlap between the %s and %s: %d" %(name_db1, name_db2, num_overlap)


def get_unique_variants(filename):
    varFile = open(filename, "r")
    variants = []
    line_num = 0
    for line in varFile:
        line_num += 1
        this_variant = line.strip().split("\t")[:4]
        alt = this_variant[-1]
        for each_alt in alt.split(","):
            branch_variant = this_variant[0:3] + [each_alt]
            if not check_variant_exist(branch_variant, variants):
                variants.append(branch_variant)
    return variants

def check_variant_exist(v, variants):
    for variant in variants:
        if v == variant:    
            return True
        elif variant_equal(v, variant):
            return True
    return False

def variant_equal(v1, v2):
    " return (edited1, edited2) "
    chr1, pos1, ref1, alt1 = v1
    chr2, pos2, ref2, alt2 = v2
    pos1 = int(pos1)
    pos2 = int(pos2)
    if chr1 != chr2:
        return False

    if (len(ref1) - len(alt1)) != (len(ref2) - len(alt2)):
        return False

    if len(ref2)>100 or len(ref1)>100:
        return False

    # make sure that v1 is upstream of v2
    if pos1 > pos2:
        return variant_equal(v2, v1)

    # lift coordinates and make everything 0-based
    if chr1 == "13":
        seq = BRCA2
        pos1 = pos1 -1 - BRCA2_START
        pos2 = pos2 -1 - BRCA2_START
    elif chr1 == "17":
        seq = BRCA1
        pos1 = pos1 - 1 - BRCA1_START
        pos2 = pos2 - 1 - BRCA1_START
    else:
        assert(False)

    # replace vcf ref string with alt string
    edited_v1 = seq[0:pos1]+alt1+seq[pos1+len(ref1):]
    edited_v2 = seq[0:pos2]+alt2+seq[pos2+len(ref2):]
    return edited_v1 == edited_v2

if __name__ == "__main__":
    main()
