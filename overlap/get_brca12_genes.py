import pysam
"""
the brca position come from UCSC genome browser and is 0-based start, 1-based end
brca1 location: chr17:41,196,312-41,277,500 81189bp
brca2 location: chr13:32,889,617-32,973,809 84193bp

this file only needs to be used once, to generate the brca1.txt and brca2.txt which is a txt file
of region covering either brca1 or brca2 genes.

"""

PATH = "/hive/groups/cgl/brca/phase1/data"
REF_PATH = PATH + "/brca_references"

def main():
    extract_brca_region_sequence()

def extract_brca_region_sequence():
    chr13 = pysam.Fastafile(REF_PATH + "/chr13.fa")
    chr17 = pysam.Fastafile(REF_PATH + "/chr17.fa")
    f13 = open(REF_PATH + "/chr13.txt", "w")
    f17 = open(REF_PATH + "/chr17.txt", "w")
    f13.write(chr13.fetch(reference="chr13", start=0, end=150000000))
    f17.write(chr17.fetch(reference="chr17", start=0, end=100000000))
    f13.close()
    f17.close()
    f13 = open(REF_PATH + "/chr13.txt", "r")
    f17 = open(REF_PATH + "/chr17.txt", "r")
    chr13 = f13.read()
    chr17 = f17.read()
    brca2 = open("brca2.txt", "w")
    brca1 = open("brca1.txt", "w")
    brca1.write(chr17[41100000:41300000])
    brca2.write(chr13[32800000:33000000])

if __name__ == "__main__":
    main()
    