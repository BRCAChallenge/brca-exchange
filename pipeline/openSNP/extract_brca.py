"""
extract SNPs from BRCA regions (including all genome assembly versions)

BRCA1: 
chr17:43,045,629-43,125,483 (GRCh38/hg38)
chr17:41,196,311-41,196,330 (GRCh37/hg19)
chr17:38,449,840-38,530,994 (NCBI36/hg18)
chr17:38,449,844-38,530,657 (NCBI35/hg17)
chr17:41,570,860-41,650,551 (NCBI34/hg16)

BRCA2:
chr13:32,315,474-32,400,266 (GRCh38/hg38)
chr13:32,889,617-32,973,809 (GRCh37/hg19)
chr13:31,787,617-31,871,809 (NCBI36/hg18)
chr13:31,787,617-31,871,805 (NCBI35/hg17)
chr13:30,687,617-30,771,805 (NCBI34/hg16)
"""

import glob



SOURCE = ["23andme", "ancestry", "ftdna"]
POSITION = {"build38": {"BRCA1": [17, 43045629, 43125483],
                        "BRCA2": [13, 32315474, 32400266]}, 
            "build37": {"BRCA1": [17, 41196311, 41196330],
                        "BRCA2": [13, 32889617, 32973809]}, 
            "build36": {"BRCA1": [17, 38449840, 38530994],
                        "BRCA2": [13, 31787617, 31871809]}, 
            "build35": {"BRCA1": [17, 38449844, 38530657],
                        "BRCA2": [13, 31787617, 31871805]}, 
            "build34": {"BRCA1": [17, 41570860, 41650551],
                        "BRCA2": [13, 30687617, 30771805]}}



def main():
    extract_23andme_brca()


def extract_23andme_brca():
    """"extract BRCA region from 23andme files by:
        1. check assembly version, if not 38, convert to 38
        2. discard all positions outside of BRCA region"""
    folder = "data/organized_opensnp_filedump/23andme/"
    for filename in glob.glob(folder + "*"):
        print filename
        f_in = open(filename, "r")
        lines = f_in.readlines()
        if not lines[0].startswith("#"):
            print "binary file"
            continue
        
        # write out header of new file
        header = [line for line in lines if line.startswith("#")]
        source = header[0]
        column_names = header[-1]
        header = "".join(header)
        build_index = header.index("human assembly build")
        version = "# " + header[build_index: build_index + 23]
        f_out = open("temp", "w")
        f_out.write(source)
        f_out.write(version + "\n")
        f_out.write(column_names)
        

        # write out body of new file with only BRCA region converted to GRCh38





if __name__ == "__main__":
    main()
