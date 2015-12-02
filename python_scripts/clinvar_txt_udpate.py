"""
This script updates the current BRCA csv file with the newest clinvar txt file
1. add extra columns to clinvar txt file: GenomeCoor, HGVS_cDNA, HGVS_protein

"""
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pygr.seqdb import SequenceFileDB

CLINVAR1 = "../../../ClinVar_XML_Parsing/ClinVar_BRCA_Nov_2015.txt"
CLINVAR2 = "../../../ClinVar_XML_Parsing/ClinVar_BRCA_Nov_2015_added_columns.txt"
CLINVAR3 = "../../../ClinVar_XML_Parsing/ClinVar_BRCA_Nov_2015_added_columns_merged_repeats.txt"

CURRENT_CSV = "../data/merge/merged_v4.csv"
UPDATED_CSV = "../data/merge/merged_v5.csv"
UPDATED_CSV2 = "../data/merge/merged_v6.csv"

def main():
    #add_Genome_Coor()
    variant_dict = remove_repeats()
    update_current_csv_with_new_variants(variant_dict)
    # update_current_csv_with_new_variant_annotation(variant_dict)

def remove_repeats():
    "remove duplicates by collapsing the Origin and ClinicalSignificance columns "
    f_in = open(CLINVAR2, "r")
    f_out = open(CLINVAR3, "w")
    variant_dict = {}
    num_repeats = 0
    num_conflicting_repeats = 0
    for line in f_in:
        if "#" in line:
            f_out.write(line)
        else:
            items = line.strip().split("\t")
            Genome_Coor = items[29]
            new_origins = set(items[11].split(";"))
            new_pathos = set(items[5].split(";"))
            if Genome_Coor in variant_dict.keys():
                num_repeats += 1
                existing_origins = variant_dict[Genome_Coor][0]
                existing_pathos = variant_dict[Genome_Coor][1]
                if (not new_origins <= existing_origins) or (not new_pathos <= existing_pathos):
                    num_conflicting_repeats += 1
                    if new_origins <= existing_origins:
                        pass
                    else:
                        variant_dict[Genome_Coor][0] = new_origins | existing_origins
                    if new_pathos <= existing_pathos:
                        pass
                    else:
                        variant_dict[Genome_Coor][1] = new_pathos | existing_pathos
            else:
                variant_dict[Genome_Coor] = [new_origins, new_pathos]
    print "number of repeats: ", num_repeats
    print "number of repeats with conflicting result", num_conflicting_repeats
    f_in.close()
    f_out.close()
    return variant_dict

def update_current_csv_with_new_variant_annotation(d):
    f_in = open(CURRENT_CSV, "r")
    f_out = open(UPDATED_CSV, "w")
    newly_identified_variants = 0
    for line in f_in:

        if "Gene_symbol" in line:
            f_out.write(line)
        else:
            items = line.strip().split(",")
            Source = items[0]
            Genome_Coor = items[2]
            if Genome_Coor in d.keys():
                items[28] = ";".join(list(d[Genome_Coor][0]))
                items[29] = ";".join(list(d[Genome_Coor][1]))
                if "ClinVar" not in Source:
                    newly_identified_variants += 1
                    items[0] = Source + "|ClinVar"
            else:
                items[28] = items[29] = "-"
            f_out.write(",".join(items) + "\n")
    print "variants that previously are not in clinvar but is now in clinvar: ", newly_identified_variants
    f_in.close()
    f_out.close()

def update_current_csv_with_new_variants(new_variants):
    f_in = open(UPDATED_CSV, "r")
    f_out = open("../data/merge/temp", "w")
    Genome_Coors = set()
    for line in f_in:
        if "Gene_symbol" not in line:
            Genome_Coor = line.strip().split(",")[2]
            Genome_Coors.add(Genome_Coor)
    num_new_variants = 0
    for key, value in new_variants.iteritems():
        if key in Genome_Coors:
            pass
        else:
            num_new_variants += 1
            items = ["-"] * 61
            items[0] = "ClinVar"
            items[2] = key
            items[28] = ";".join(list(value[0]))
            items[29] = ";".join(list(value[1]))
            f_out.write(",".join(items) + "\n")
    print "novel variants: ", num_new_variants

def add_Genome_Coor():
    f_in = open(CLINVAR1, "r")
    f_out = open(CLINVAR2, "w")

    for line in f_in:
        if "#" in line:
            items = line.strip().split("\t")
            items += ["GenomeCoordinate", "HGVS_cDNA", "HGVS_protein"]
            f_out.write(",".join(items) + "\n")
        else:
            items = line.strip().split("\t")
            HGVS = items[2]
            if "NM_" in HGVS: # ignore the non-HGVS nomenclature which is about 200/16000 cases
                # remove the (BRCA1) string so that hgvs module can parse the HGVS string
                if "BRCA" in HGVS:
                    stuff = HGVS.split(" ")[0].split(":")
                    HGVS_c = stuff[0][0:-7] + ":" + stuff[1]
                else:
                    HGVS_c = HGVS

                try: # only get HGVS_p from txt file if it exists
                    HGVS_p = HGVS.split(" ")[1][1:-1]
                except:
                    HGVS_p = ""

                try:
                    GenomeCoor = HGVS_to_GenomeCoor(HGVS_c)
                except:
                    GenomeCoor = None

                if GenomeCoor is not None:
                    items += [GenomeCoor, HGVS_c, HGVS_p]
                    f_out.write("\t".join(items) + "\n")
    f_in.close()
    f_out.close()



def HGVS_to_GenomeCoor(HGVS):
    """use counsyl pyhgvs for this"""
    genome = SequenceFileDB('../data/hg19.fa')
    refGene = "../data/BRCA12.refGene.txt"
    with open(refGene) as infile:
        transcripts = pyhgvs_utils.read_transcripts(infile)
    def get_transcript(name):
        return transcripts.get(name)
    chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(
        HGVS, genome, get_transcript=get_transcript)
    genome_coordinate = chrom + ":" + str(offset) + ":" + ref + ">" + alt
    return genome_coordinate



if __name__ == '__main__':
    main()

