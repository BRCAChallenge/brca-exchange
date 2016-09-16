"""extract BRCA region from 23andme and ancestryDNA files by:
   1. check assembly version, if not 38, convert positions to 38
   2. discard all positions outside of BRCA region
    Note: about 2000 files from 23andme, 200 files from ancestryDNA, 200 files from ftdna
    ftdna files are not processed because its human reference build version is not specified
"""

import glob
import pyliftover
import pdb
import os


SOURCE = ["23andme", "ancestry", "ftdna"]

# preload all the pyliftover functions because the function needs to download from internet
LIFT_MAP = {"37": pyliftover.LiftOver('hg19', 'hg38'),
               "36": pyliftover.LiftOver('hg18', 'hg38'),
               # pyliftover refuses to translate from hg17 to hg38
               # therefore it's done in two steps hg17 -> hg19 -> hg38
               "35": [pyliftover.LiftOver('hg17', 'hg19'), pyliftover.LiftOver('hg19','hg38')],
               "34": pyliftover.LiftOver('hg16', 'hg38')}

BRCA_BOUNDARY = {"38": {"chr17": [43045629, 43125483],
                        "chr13": [32315474, 32400266]}, 
                 "37": {"chr17": [41196312, 41277500],
                        "chr13": [32889617, 32973809]}, 
                 "36": {"chr17": [38449840, 38530994],
                        "chr13": [31787617, 31871809]}, 
                 "35": {"chr17": [38449844, 38530657],
                        "chr13": [31787617, 31871805]}, 
                 "34": {"chr17": [41570860, 41650551],
                        "chr13": [30687617, 30771805]}}


def genomic_coordinate_update(version, chrm, position):
    if version == "35":
        lift_tool1 = LIFT_MAP[version][0]
        lift_tool2 = LIFT_MAP[version][1]
        intermediate_position = lift_tool1.convert_coordinate(chrm, position)[0][1]
        return lift_tool2.convert_coordinate(chrm, position)[0][1]
    else:
        lift_tool = LIFT_MAP[version]
        return lift_tool.convert_coordinate(chrm, position)[0][1]


def main():
    outpath = "data/brca_openSNP/"
    try: 
        os.mkdir(outpath)
    except OSError:
        pass
    extract_23andme_brca(outpath)
    extract_ancestry_brca(outpath)


def extract_ancestry_brca(outpath):
    folder = "data/organized_opensnp_filedump/ancestry/"
    for index, filename in enumerate(glob.glob(folder + "*")):
        filename_short = filename.split("/")[-1]
        print index, filename_short
        f_in = open(filename, "r")
        lines = f_in.readlines()
        if not lines[0].startswith("#"):
            print "skipping binary file"
            continue
        
        # write out header of new file
        header = [line for line in lines if line.startswith("#")]
        column_names = header[-1]
        header = "".join(header)
        build_index = header.index("human reference build")
        version = header[build_index: build_index + 24][-2:]
        f_out = open(outpath + filename_short, "w")
        f_out.write("# ancestry\n")
        f_out.write("# human reference build 38\n")
        
        # write out body of new file with only BRCA region converted to GRCh38
        body = [line for line in lines if not line.startswith("#")]
        for line_index, line in enumerate(body):
            items = line.strip().split("\t")
            if line_index == 0:
                items.pop(-1)
                items[-1] = "Alleles"
                f_out.write("#" + "\t".join(items) + "\n")
                continue
            chrm = "chr" + items[1]
            position = int(items[2])
            if chrm == "chr17" or chrm == "chr13":
                [brca_start, brca_end] = BRCA_BOUNDARY[version][chrm]
                if position >= brca_start and position <= brca_end:
                    if version == "38":
                        f_out.write(line)
                    else:
                        new_pos = int(genomic_coordinate_update(version, chrm, position))
                        items[2] = str(new_pos)
                        items[-2] = items[-2] + items[-1]
                        items.pop(-1)
                        f_out.write("\t".join(items) + "\n")     
                else:
                    continue
            else:
                continue
        f_in.close()
        f_out.close()


def extract_23andme_brca(outpath):
    
    folder = "data/organized_opensnp_filedump/23andme/"
    for index, filename in enumerate(glob.glob(folder + "*")):
        filename_short = filename.split("/")[-1]
        print index, filename_short
        f_in = open(filename, "r")
        lines = f_in.readlines()
        if not lines[0].startswith("#"):
            print "skipping binary file"
            continue
        
        # write out header of new file
        header = [line for line in lines if line.startswith("#")]
        column_names = header[-1]
        header = "".join(header)
        build_index = header.index("human assembly build")
        version = header[build_index: build_index + 23][-2:]
        f_out = open(outpath + filename_short, "w")
        f_out.write("# 23andme\n")
        f_out.write("# human assembly build 38\n")
        f_out.write(column_names)
        
        # write out body of new file with only BRCA region converted to GRCh38
        body = [line for line in lines if not line.startswith("#")]
        for line in body:
            items = line.strip().split("\t")
            chrm = "chr" + items[1]
            position = int(items[2])
            if chrm == "chr17" or chrm == "chr13":
                [brca_start, brca_end] = BRCA_BOUNDARY[version][chrm]
                if position >= brca_start and position <= brca_end:
                    if version == "38":
                        f_out.write(line)
                    else:
                        new_pos = int(genomic_coordinate_update(version, chrm, position))
                        items[2] = str(new_pos)
                        f_out.write("\t".join(items) + "\n")     
                else:
                    continue
            else:
                continue
        f_in.close()
        f_out.close()


if __name__ == "__main__":
    main()
