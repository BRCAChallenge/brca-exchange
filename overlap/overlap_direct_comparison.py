import glob


# remember the last "/", the PATH won't work without the trailing "/"`
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
        db_name = file.split("/")[-1].split(".")[0]
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
