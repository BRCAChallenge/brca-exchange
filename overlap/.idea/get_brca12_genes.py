import pysam
"""
NOTE: the position come from UCSC genome browser and is 0-based start, 1-based end
brca1 location: chr17:41,196,312-41,277,500 81189bp
brca2 location: chr13:32,889,617-32,973,809 84193bp
"""


def main():
    # BRCA2_START = 32800000
    # BRCA1_START = 41100000
    #
    # chr13 = open("data/brca2.txt", "r")
    # chr13 = chr13.read()
    # print chr13[32890539-BRCA2_START-1]
    # print chr13[32890543-BRCA2_START-1:32890543-BRCA2_START-1+4]
    #
    # chr17 = open("data/brca1.txt", "r")
    # chr17 = chr17.read()
    # print chr17[41276033-BRCA1_START-1:41276033-BRCA1_START-1+10]
    extract_brca_region_sequence()



def turn_fasta_to_txt(filename, chrm, start, end):
    f = pysam.Fastafile(filename)
    text = open("data/{0}.txt".format(chrm), "w")
    text.write(f.fetch(reference=chrm, start=start, end=end))

def extract_brca_region_sequence():
    chr13 = pysam.Fastafile("references/chr13.fa")
    chr17 = pysam.Fastafile("references/chr17.fa")
    f13 = open("references/chr13.txt", "w")
    f17 = open("references/chr17.txt", "w")
    f13.write(chr13.fetch(reference="chr13", start=0, end=150000000))
    f17.write(chr17.fetch(reference="chr17", start=0, end=100000000))
    f13.close()
    f17.close()
    f13 = open("references/chr13.txt", "r")
    f17 = open("references/chr17.txt", "r")
    chr13 = f13.read()
    chr17 = f17.read()
    brca2 = open("references/brca2.txt", "w")
    brca1 = open("references/brca1.txt", "w")
    brca1.write(chr17[41100000:41300000])
    brca2.write(chr13[32800000:33000000])


if __name__ == "__main__":
    main()
    