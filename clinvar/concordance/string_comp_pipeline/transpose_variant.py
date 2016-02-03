"""
this script turns data into a variant centric view
each variant (as by uniq_id) occupies one row, and all the submitters
who reported this variant as well the reported pathogenicity are shown in this row
"""
import pandas as pd
import numpy as np
from pprint import pprint
import string_comp

SUBMITTERS = ["Invitae", "GeneDx", "Ambry", "Emory", "Children",
              "SCRP", "Counsyl"]

COLUMNS = ["uniq_id", "HGVS"] + SUBMITTERS

IN_FILE = "../BRCA_selectedLabs_only/BRCA.pre-processed.add_uniq_id"
IN_FILE1 = IN_FILE + ".add_concordance"
OUT_FILE = "../BRCA_selectedLabs_only/BRCA.transposed_by_uniq_id"

def main():
    transpose()
    add_concordance()

def transpose():
    df = pd.read_csv(IN_FILE, sep="\t")
    uniq_id_dict = {}

    for index, row in df.iterrows():
        this_submission = "|".join(
                [row.ClinicalSignificance, row.DateLastUpdated, row.DateCreated, row.SCV])
        if row.uniq_id not in uniq_id_dict.keys():
            uniq_id_dict[row.uniq_id] = {"HGVS": set([row.HGVS]),
                                         "Invitae": [],
                                         "GeneDx": [],
                                         "Ambry": [],
                                         "Emory": [],
                                         "Children": [],
                                         "SCRP": [],
                                         "Counsyl": []}
            for submitter_abrev in SUBMITTERS:
                if submitter_abrev.upper() in row.Submitter.upper():
                    uniq_id_dict[row.uniq_id][submitter_abrev].append(this_submission)
        else:
            uniq_id_dict[row.uniq_id]["HGVS"].add(row.HGVS)
            for submitter_abrev in SUBMITTERS:
                if submitter_abrev.upper() in row.Submitter.upper():
                    uniq_id_dict[row.uniq_id][submitter_abrev].append(this_submission)


    # transform the format to fit pandas.dataframe
    row_list = []
    for uniq_id, value in uniq_id_dict.iteritems():
        this_dict = {"uniq_id": uniq_id}
        this_dict["HGVS"] = "|".join(list(value["HGVS"]))
        for submitter_abrev in SUBMITTERS:
            this_dict[submitter_abrev] = ",".join(value[submitter_abrev])
        row_list.append(this_dict)

    transposed_df = pd.DataFrame(row_list)
    transposed_df = transposed_df[COLUMNS]
    transposed_df.to_csv(OUT_FILE, sep="\t", dtype=str, index=False)

def add_concordance():
    df = pd.read_csv(OUT_FILE, sep="\t", dtype=str)
    n_submitters = []
    n_submissions = []
    concordance = []
    for index, row in df.iterrows():
        n_submission = 0
        n_submitter = 0
        reported_results = set()
        for submitter in SUBMITTERS:
            if not pd.isnull(row[submitter]):
                n_submitter += 1
                for report in row[submitter].split(","):
                    reported_results.add(report.split("|")[0])
                n_submission += len(row[submitter].split(","))
        n_submitters.append(n_submitter)
        n_submissions.append(n_submission)
        concordance.append(decide_concordance(reported_results))
    df["num_submitters"] = pd.DataFrame(n_submitters)
    df['concordance'] = pd.DataFrame(concordance)
    df.to_csv(OUT_FILE+".add_concordance", sep="\t", index=False, dtype=str)


def decide_concordance(patho_set):
    if len(patho_set) == 1:
        return True
    else:
        concordant = True
        patho_list = list(patho_set)
        while len(patho_list) > 0:
            first_item = patho_list.pop(0)
            for each_of_rest in patho_list:
                if is_discordant(first_item, each_of_rest):
                    concordant = False             
        return concordant

def is_discordant(term1, term2):
    """assess concordinance based on invitae poster criteria"""
    if "PATHO" in (term1.upper() + term2.upper()):
        if ("PATHO" in term1.upper()) and ("PATHO" in term2.upper()):
            return False
        else:
            return True
    else:
        return False



if __name__ == "__main__":
    main()


