def main():
    f = open("../data/allVcf/bic_brca12.sorted.vcf", "r")
    for line in f:
        items = line.strip().split("\t| +")
        try:
            x = items[8]
            print line
            print items
        except IndexError:
            pass


if __name__ == "__main__":
    main()

