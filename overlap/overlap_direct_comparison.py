import glob



"""
    result:
unique variants in ClinVar:6065
overlap between ClinVar and ex_LOVD: 231
overlap between ClinVar and LOVD: 1558
overlap between ClinVar and 1000_Genomes: 341
overlap between ClinVar and ExAC: 1191
overlap between ClinVar and BIC: 3203
overlap between ClinVar and UMD: 1395

unique variants in ex_LOVD:288
overlap between ex_LOVD and ClinVar: 231
overlap between ex_LOVD and LOVD: 265
overlap between ex_LOVD and 1000_Genomes: 83
overlap between ex_LOVD and ExAC: 177
overlap between ex_LOVD and BIC: 259
overlap between ex_LOVD and UMD: 209

unique variants in LOVD:3129
overlap between LOVD and ClinVar: 1558
overlap between LOVD and ex_LOVD: 265
overlap between LOVD and 1000_Genomes: 352
overlap between LOVD and ExAC: 911
overlap between LOVD and BIC: 1479
overlap between LOVD and UMD: 1168

unique variants in 1000_Genomes:4351
overlap between 1000_Genomes and ClinVar: 341
overlap between 1000_Genomes and ex_LOVD: 83
overlap between 1000_Genomes and LOVD: 352
overlap between 1000_Genomes and ExAC: 520
overlap between 1000_Genomes and BIC: 317
overlap between 1000_Genomes and UMD: 363

unique variants in ExAC:3738
overlap between ExAC and ClinVar: 1191
overlap between ExAC and ex_LOVD: 177
overlap between ExAC and LOVD: 911
overlap between ExAC and 1000_Genomes: 520
overlap between ExAC and BIC: 924
overlap between ExAC and UMD: 859

unique variants in BIC:3633
overlap between BIC and ClinVar: 3203
overlap between BIC and ex_LOVD: 259
overlap between BIC and LOVD: 1479
overlap between BIC and 1000_Genomes: 317
overlap between BIC and ExAC: 924
overlap between BIC and UMD: 1178

unique variants in UMD:3854
overlap between UMD and ClinVar: 1395
overlap between UMD and ex_LOVD: 209
overlap between UMD and LOVD: 1168
overlap between UMD and 1000_Genomes: 363
overlap between UMD and ExAC: 859
overlap between UMD and BIC: 1178
"""


PATH = "/hive/groups/cgl/brca/phase1/data/cutoff_vcf/"


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
        db_name = file.split(".")[0].split("/")[1]
        db[db_name] = get_variant_list(file)
    return db

def print_variant_size(database, database_name):
    print "unique variants in %s:%d" %(database_name, len(set(database)))

def get_overlap(name1, name2, database):
    print "overlap between %s and %s: %d" %(name1, name2, len(set(database[name1]).intersection(set(database[name2]))))

def get_variant_list(filename):
    varFile = open(filename, "r")
    variants = []
    for line in varFile:
        items = line.split("\t")
        chr = items[0]
        pos = items[1]
        ref = items[2]
        alts = items[3].split(",")
        for alt in alts:
            this_variant = "{0}.{1}.{2}.{3}".format(chr, pos, ref, alt)
            variants.append(this_variant)
    return variants

if __name__ == "__main__":
    main()
