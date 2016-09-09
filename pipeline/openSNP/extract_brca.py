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
chr17:41,570,860-41,650,551 (NCBI34/hg16)
"""

import glob




RAW_DATA = "data/datadump_opensnp_9_8_2016/"



def main():
    file_stats()




def file_stats():
    file_dict = {"picture": [],
                 "fitbit": [],
                 "23andme": [],
                 "ancestry": [],
                 "ftdna" :[],
                 "IYG": [],
                 "decodeme": []}
    subfix = {}
    for index, filename in enumerate(glob.glob(RAW_DATA + "*")):
        
        # find out file subfix distribution
        filetype = filename.split(".")[-1]
        if filetype in subfix.keys():
            subfix[filetype] += 1
        else:
            subfix[filetype] = 1

        # seperate files based on content
        found = False
        for key in file_dict.keys():
            if key in filename:
                found = True
                file_dict[key].append(filename)
        
        if not found:
            print filename    


    print subfix
    sum = 0
    for key, value in file_dict.iteritems():
        print key
        print len(value)
        sum += len(value) 
    print sum




if __name__ == "__main__":
    main()
