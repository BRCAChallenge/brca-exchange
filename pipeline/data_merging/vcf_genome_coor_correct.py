import vcf
import glob
import string_comp

PATH = "/hive/groups/cgl/brca/phase1/data/allVcfs/"




def main():
    for file in glob.glob(PATH + "*"):
        f = open(file, "r")
        vcf_reader = vcf.Reader(f, strict_whitespace=True)
        n_wrong = 0
        n_total = 0
        for record in vcf_reader:
            n_total += 1
            v = [record.CHROM, record.POS, record.REF, record.ALT]
            if not string_comp.ref_correct(v, version="hg19"):
                n_wrong += 1
        print "in {0}\nwrong: {1}, total: {2}".format(
            file, n_wrong, n_total) 








if __name__ == "__main__":
    main()
