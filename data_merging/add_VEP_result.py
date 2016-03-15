#!/usr/bin/env python
"""add VEP result to merged data"""
import argparse
import pandas as pd
import json
import pickle

VEP_OUTPUT = "/cluster/home/mollyzhang/release1.0/data/VEP/vep_output_3_3_2016.vcf"
VEP_FIELDS = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene',
              'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON',
              'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position',
              'Protein_position', 'Amino_acids', 'Codons',
              'Existing_variation', 'DISTANCE', 'STRAND', 'SYMBOL_SOURCE',
              'HGNC_ID', 'TSL', 'APPRIS', 'REFSEQ_MATCH', 'SIFT', 'PolyPhen',
              'GMAF', 'AFR_MAF', 'AMR_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF',
              'AA_MAF', 'EA_MAF', 'ExAC_MAF', 'ExAC_Adj_MAF', 'ExAC_AFR_MAF',
              'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF', 'ExAC_NFE_MAF',
              'ExAC_OTH_MAF', 'ExAC_SAS_MAF', 'CLIN_SIG', 'SOMATIC', 'PHENO',
              'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS',
              'MOTIF_SCORE_CHANGE']
#UNWANTED = ["Allele", "SYMBOL", "Gene", "Feature_type", "Feature", "HGVSc",
#            "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
#            "SYMBOL_SOURCE", "HGNC_ID", "GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF##,
#            "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF", "ExAC_MAF", 'ExAC_Adj_MA##',
#            'ExAC_AFR_MAF', 'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF',
#            'ExAC_NFE_MAF', 'ExAC_OTH_MAF', 'ExAC_SAS_MAF', "CLIN_SIG"]
UNWANTED = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene',
            'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON',
            'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position',
            'Protein_position', 'Amino_acids', 'Codons',
            'Existing_variation', 'DISTANCE', 'STRAND', 'SYMBOL_SOURCE',
            'HGNC_ID', 'TSL', 'APPRIS', 'REFSEQ_MATCH',
            'GMAF', 'AFR_MAF', 'AMR_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF',
            'AA_MAF', 'EA_MAF', 'ExAC_MAF', 'ExAC_Adj_MAF', 'ExAC_AFR_MAF',
            'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF', 'ExAC_NFE_MAF',
            'ExAC_OTH_MAF', 'ExAC_SAS_MAF', 'CLIN_SIG', 'SOMATIC', 'PHENO',
            'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS',
            'MOTIF_SCORE_CHANGE']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        default="/hive/groups/cgl/brca/release1.0/merged.csv")
    parser.add_argument("-o", "--output",
                        default="/hive/groups/cgl/brca/release1.0/merged_withVEP_cleaned.csv")
    args = parser.parse_args()
    vep_result_dict = save_VEP_to_dict()
    temp_dump = open("temp_dump", "w")
    temp_dump.write(pickle.dumps(vep_result_dict))
    temp_dump.close()
    #vep_result_dict = pickle.loads(open("temp_dump", "r").read())
    write_to_file(vep_result_dict, args.input, args.output)

def write_to_file(vep_result_dict, inputFile, outputFile):
    print "write to file {0}".format(outputFile)
    f_in = open(inputFile, "r")
    f_out = open(outputFile, "w")
    line_num = 0
    for field in UNWANTED:
        VEP_FIELDS.remove(field)
    for line in f_in:
        line_num += 1
        if line_num %1000 == 0:
            print line_num
        items = line.strip().split(",")
        if line_num == 1:
            vep_fields = [i+"(VEP)" for i in VEP_FIELDS]
            items += vep_fields
        else:
            genome_coor = items[2]
            additional_items = []
            coordinate_this_variant = ""
            for coordinate in genome_coor.split("|"):
                if coordinate in vep_result_dict.keys():
                    coordinate_this_variant = coordinate
            if coordinate_this_variant == "":
                for column in VEP_FIELDS:
                    additional_items.append("-")
            else:
                for column in VEP_FIELDS:
                    this_cell = vep_result_dict[coordinate_this_variant][column].replace(", ","|")
                    additional_items.append(this_cell)
            items += additional_items
        new_line = ",".join(items) + "\n"
        f_out.write(new_line)



def save_VEP_to_dict():
    print("processing VEP output")
    f = open(VEP_OUTPUT, "r")
    vep_dict = {}
    line_no = 0
    for line in f:
        if "#" in line:
            continue
        line_no += 1
        if line_no%10 == 0: 
            print(line_no)
        items = line.replace(".", "-").strip().split("\t")
        genome_coor = "chr{0}:{1}:{2}>{3}".format(items[0],
                                                  items[1],
                                                  items[3],
                                                  items[4])
        if genome_coor not in vep_dict:
            if items[-1] != "-":
                vep_dict[genome_coor] = majority_vote_VEP(items[-1])
        else:
            raise Exception("bad time")
    return vep_dict

def collapse_VEP_info(info):
    """collapse various VEP results into one by giving each one a weight"""
    infos = info[4:].split(",")
    row_list = []
    for each_info in infos:
        values = each_info.split("|")
        result_dict = dict(zip(VEP_FIELDS, values))
        row_list.append(result_dict)
    df = pd.DataFrame(row_list)
    collapsed = {}
    for field in df.columns:
        this_field = df[field].value_counts()
        if len(this_field) > 1:
            value_list = []
            for key, value in df[field].value_counts().iteritems():
                value_list.append([key, "{0}/{1}".format(value, str(len(infos)))])
            collapsed[field] = json.dumps(value_list)
        else:
            collapsed[field] = this_field.index[0]
    return collapsed

def majority_vote_VEP(info):
    """take the majority vote VEP"""
    infos = info[4:].split(",")
    row_list = []
    for each_info in infos:
        values = each_info.split("|")
        result_dict = dict(zip(VEP_FIELDS, values))
        for key in UNWANTED:
            del result_dict[key]
        row_list.append(result_dict)
    df = pd.DataFrame(row_list)
    majority = {}
    for field in df.columns:
        majority[field] = df[field].value_counts().index[0]
    return majority



if __name__ == "__main__":
    main()

