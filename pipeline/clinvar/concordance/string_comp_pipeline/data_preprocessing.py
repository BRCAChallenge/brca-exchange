"""
Data preprocessing:
1. save a log of information about original data file to filename.log
2. remove repetitive rows from invitae with same HGVS string (save the older submission)
2. for row that do have none/wrong genomic coordinate, translate genome coor from HGVS

"""

import operator
import pandas as pd
import string_comp
import pickle
import collections

import HGVS_conversion


FILENAME = "../BRCA_selectedLabs_only/ClinVarBRCA.selectedLabsOnly.txt"
FILENAME1 = "../BRCA_selectedLabs_only/ClinVarBRCA.selectedLabsOnly.remove_invitae_dup.txt"

DF = pd.read_csv(FILENAME1, sep="\t", dtype=str)
NEW_FILE = "../BRCA_selectedLabs_only/BRCA.pre-processed"
ERROR = "../BRCA_selectedLabs_only/BRCA.wrong_genome_Coor"



def main():
    remove_invitae_dup()
    save_log(FILENAME1)
    add_genome_coor()

def remove_invitae_dup():
    df = pd.read_csv(FILENAME, sep="\t", dtype=str)
    df_no_invitae = df.loc[df.Submitter != "Invitae", :]
    df_invitae = df.loc[df.Submitter == "Invitae", :]

    duplicate_HGVS = [i for i, count in collections.Counter(df_invitae.HGVS).items() if count > 1]
    single_HGVS = [i for i, count in collections.Counter(df_invitae.HGVS).items() if count == 1]


    row_list_dup = []
    for dup in duplicate_HGVS:
        for index, row in df_invitae[df_invitae.HGVS == dup].iterrows():
            if row.SCV != "SCV000000000":
                row_list_dup.append(row.to_dict())

    row_list_single=[]
    for single in single_HGVS:
        for index, row in df_invitae[df_invitae.HGVS == single].iterrows():
            row_list_single.append(row.to_dict())

    df_invitae_de_dup = pd.DataFrame(row_list_dup)
    df_invitae_de_dup = df_invitae_de_dup[df.columns]
    df_invitae_single = pd.DataFrame(row_list_single)
    df_invitae_single = df_invitae_single[df.columns]

    frames = [df_no_invitae, df_invitae_single, df_invitae_de_dup]
    new_df = pd.concat(frames)
    new_df.to_csv(FILENAME1, sep="\t", dtype=str, index=False)


def save_log(filename):
    log = open(filename + ".log", "w")
    log.write(filename + "\n")
    log.write("number of rows: {0}, \n".format(len(DF)))
    log.write("number of repetitive rows: {0}\n".format(
        len(DF) - len(DF.drop_duplicates())))
    log.write("number of unique variant ids: {0}\n".format(len(set(DF.VariantID))))
    log.write("number of unique HGSV strings: {0}\n\n".format(len(set(DF.HGVS))))

    write_submitter_and_clinsigs(log, DF)

def add_genome_coor():
    temp_genome_coor = (DF['Chrom'] + "_" + DF['Pos'] + "_" +
                     DF['Ref'].replace("-", "") + "_" + DF['Alt'].replace("-", ""))
    columns = DF.columns
    error_row_list = []
    genome_coor = []
    for index, row in DF.iterrows():
        variant = temp_genome_coor[index].split("_")
        if (("None" in variant) or (not string_comp.ref_correct(variant))):
            genome_coor.append(HGVS_conversion.HGVS_to_genome_coor(row.HGVS))
        else:
            genome_coor.append(temp_genome_coor[index])
        if "None" not in variant and not string_comp.ref_correct(variant):
            error_row_list.append(row.to_dict())

    DF['genome_coor'] = pd.Series(genome_coor)
    DF.to_csv(NEW_FILE, sep="\t", index=False)
    df_error = pd.DataFrame(error_row_list)
    df_error = df_error[columns]
    df_error.to_csv(ERROR, sep="\t", index=False, dtype=str)


def write_submitter_and_clinsigs(log, df):
    submitters = {}
    clinsigs = {}
    for index, row in df.iterrows():
        if submitters.get(row.Submitter) is not None:
            submitters[row.Submitter] += 1
        else: submitters[row.Submitter] = 1
        if clinsigs.get(row.ClinicalSignificance) is not None:
            clinsigs[row.ClinicalSignificance] += 1
        else: clinsigs[row.ClinicalSignificance] = 1
    log.write("num_submission\tsubmitter\n")
    write_dict_in_order(log, submitters)
    log.write("case_num\tclinical_significance\n")
    write_dict_in_order(log, clinsigs)

def write_dict_in_order(file, dict):
    sorted_dict = sorted(dict.items(), key=operator.itemgetter(1))
    for key in sorted_dict:
        file.write(str(key[1]) + "\t" + key[0] + "\n")
    file.write(str(sum(dict.values())) + "\t" + "Total\n\n")




if __name__=="__main__":
    main()

