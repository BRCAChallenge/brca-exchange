import glob
import vcfFindDifferentButSame as max


"""
    improved algorithm for comparing two variants:
    if different chromosome? different variant
    if position of two variants is more than 300 bp apart? different variant
    if:
        restore the common sequence between the two variants (e.g. start with earlier variant ending with later variant)
        compare the common sequence, if different: different variant
    else: same variant

    result:
unique variants in ClinVar:6062
overlap between the ClinVar and ex_LOVD: 274
overlap between the ClinVar and LOVD: 1766
overlap between the ClinVar and 1000_Genomes: 386
overlap between the ClinVar and ExAC: 1399
overlap between the ClinVar and BIC: 3512
overlap between the ClinVar and UMD: 1562

unique variants in ex_LOVD:288
overlap between the ex_LOVD and ClinVar: 274
overlap between the ex_LOVD and LOVD: 265
overlap between the ex_LOVD and 1000_Genomes: 84
overlap between the ex_LOVD and ExAC: 198
overlap between the ex_LOVD and BIC: 259
overlap between the ex_LOVD and UMD: 210

unique variants in LOVD:3072
overlap between the LOVD and ClinVar: 1766
overlap between the LOVD and ex_LOVD: 265
overlap between the LOVD and 1000_Genomes: 354
overlap between the LOVD and ExAC: 1009
overlap between the LOVD and BIC: 1483
overlap between the LOVD and UMD: 1176

unique variants in 1000_Genomes:4351
overlap between the 1000_Genomes and ClinVar: 386
overlap between the 1000_Genomes and ex_LOVD: 84
overlap between the 1000_Genomes and LOVD: 354
overlap between the 1000_Genomes and ExAC: 559
overlap between the 1000_Genomes and BIC: 318
overlap between the 1000_Genomes and UMD: 364

unique variants in ExAC:3738
overlap between the ExAC and ClinVar: 1389
overlap between the ExAC and ex_LOVD: 197
overlap between the ExAC and LOVD: 1002
overlap between the ExAC and 1000_Genomes: 558
overlap between the ExAC and BIC: 1011
overlap between the ExAC and UMD: 945

unique variants in BIC:3630
overlap between the BIC and ClinVar: 3512
overlap between the BIC and ex_LOVD: 259
overlap between the BIC and LOVD: 1483
overlap between the BIC and 1000_Genomes: 318
overlap between the BIC and ExAC: 1017
overlap between the BIC and UMD: 1187

unique variants in UMD:3733
overlap between the UMD and ClinVar: 1561
overlap between the UMD and ex_LOVD: 211
overlap between the UMD and LOVD: 1179
overlap between the UMD and 1000_Genomes: 364
overlap between the UMD and ExAC: 950
overlap between the UMD and BIC: 1190


"""


#PATH = "/hive/groups//brca/phase1/data/cutoff_vcf/"
PATH = "/Users/Molly/PycharmProjects/brca_personal/cutoff_vcf/"
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
        print line_num, filename
        this_variant = line.strip().split("\t")[:4]
        alt = this_variant[-1]
        if "," in alt:
            alts = alt.split(",")
            for each_alt in alts:
                branch_variant = this_variant[0:3] + [each_alt]
                if not check_variant_exist(branch_variant, variants):
                    variants.append(branch_variant)
        else:
            if not check_variant_exist(this_variant, variants):
                variants.append(this_variant)
    return variants

def check_variant_exist(v, variants):
    for variant in variants:
        ### only run the comparison algorithm if two variants are different
        ### if they are the same, no need to restore the sequence ot compare
        if v == variant:
            return True
        # elif variant_pair_same(v, variant):
        #     return True
        else:
            max_result =  max.variant_seqs(v, variant)
            if max_result[0] == max_result[1]:
                return True
    return False

def variant_pair_same(v1, v2):
    chr1, pos1, ref1, alt1 = v1
    chr2, pos2, ref2, alt2 = v2
    if chr1 != chr2:
        return False
    pos1 = int(pos1)
    pos2 = int(pos2)
    # make sure that v1 is upstream of v2
    if pos1 > pos2:
        return variant_pair_same(v2, v1)

    # include the length of "alt" or "ref" of downstream variant in distance between two variants
    # v1: -------AT--------------
    # v2: ---------------CCG-----
    # the common sequence should include the section from A to G
    distance = max(pos2 - pos1 + max(len(ref2), len(alt2)),
                   max(len(ref1), len(alt1)))
    if distance > 300:
        return False

    if chr1 == "13":
        brca_string = BRCA2
        brca_start = BRCA2_START

    elif chr1 == "17":
        brca_string = BRCA1
        brca_start = BRCA1_START

    # change the coordinates to the entire chromsome to just brca1/2 genes
    pos1 = pos1 - 1 - brca_start
    pos2 = pos2 - 1 - brca_start

    ref = brca_string
    edited_v1 = ref[:pos1] + alt1 + ref[(pos1 + len(alt1)):]
    edited_v2 = ref[:pos2] + alt2 + ref[(pos2 + len(alt2)):]

    return edited_v1 == edited_v2

if __name__ == "__main__":
    main()