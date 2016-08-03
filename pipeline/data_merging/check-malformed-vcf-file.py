def main():
    f = open("../data/allVcf/no_repeats/bic.brca.no_repeats.vcf", "r")
    for line in f:
        items = line.strip().split("\t")
        try:
            x = items[8]
            print "sample column exisits"
        except IndexError:
            pass
        if len(items) != 8:
            print len(items)
            print "column numbers is not 8"
            print line
if __name__ == "__main__":
    main()

