"""add VEP result to merged data"""
import pandas as pd
import json
import pickle

MERGED_FILE = "/cluster/home/mollyzhang/release1.0/merged.csv"
OUTPUT = "/cluster/home/mollyzhang/release1.0/merged_withVEP.csv"
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



def main():
    vep_result_dict = save_VEP_to_dict()
    #temp_dump = open("temp_dump", "w")
    #temp_dump.write(pickle.dumps(vep_result_dict))
    #temp_dump.close()
    #vep_result_dict = pickle.loads(open("temp_dump", "r").read())
    write_to_file(vep_result_dict)

def write_to_file(vep_result_dict):
    f_in = open(MERGED_FILE, "r")
    f_out = open(OUTPUT, "w")
    line_num = 0
    for line in f_in:
        line_num += 1
        if line_num %100 == 0:
            print line_num
        items = line.strip().split(",")
        if line_num == 1:
            vep_fields = ["VEP_" + i for i in VEP_FIELDS]
            items += vep_fields
        else:
            genome_coor = items[2]
            if genome_coor in vep_result_dict.keys():
                additional_items = []
                for column in VEP_FIELDS:
                    this_cell = vep_result_dict[genome_coor][column].replace(", ","|")
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
                vep_dict[genome_coor] = process_VEP_info(items[-1])
        else:
            raise Exception("bad time")
    return vep_dict

def process_VEP_info(info):
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


if __name__ == "__main__":
    main()

