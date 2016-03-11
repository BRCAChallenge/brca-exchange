#!/usr/bin/env python

import argparse
import csv
import re

EMPTY = "-"
FIELDS_TO_REMOVE=["Gene_symbol(ENIGMA)", "Genomic_Coordinate(ENIGMA)",
                  "Reference_sequence(ENIGMA)", "HGVS_cDNA(ENIGMA)",
                  "BIC_Nomenclature(ENIGMA)", "Abbrev_AA_change(ENIGMA)",
                  "HGVS_protein(ENIGMA)", 
                  "HGVS(ClinVar)", "HGVS_cDNA(LOVD)", "HGVS_cDNA(exLOVD)",
                  "HGVS_protein(LOVD)", "HGVS_protein(exLOVD)",
                  "polyPhen2_result(ESP)", 
                  "BIC_Nomenclature(BIC)", "BIC_Nomenclature(exLOVD)"]
FIELDS_TO_ADD=["Gene_Symbol", "Genomic_Coordinate", "Reference_Sequence",
               "HGVS_cDNA", "BIC_Identifier", "HGVS_Protein", 
               "Protein_Change", "Allele_Frequency", 
               "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36", 
               "Source_URL", "Discordant", "Other_HGVS_cDNA"]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        default="/hive/groups/cgl/brca/release1.0/merged_withVEP_cleaned.csv")
    parser.add_argument("-o", "--output",
                        default="/hive/groups/cgl/brca/release1.0/aggregated.csv")
    args = parser.parse_args()

    csvIn = csv.DictReader(open(args.input, "r"))
    outputColumns = setOutputColumns(csvIn.fieldnames, FIELDS_TO_REMOVE,
                                     FIELDS_TO_ADD)
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter=',',
                            fieldnames=outputColumns)
    csvOut.writerow(dict((fn,fn) for fn in outputColumns))
    for row in csvIn:
        csvOut.writerow(updateRow(row, FIELDS_TO_REMOVE))

        
def setOutputColumns(fields, toRemove, toAdd):
    newFields = []
    for item in fields:
        newFields.append(item)
    for item in toRemove:
        newFields.remove(item)
    for item in toAdd:
        newFields.append(item)
    return(newFields)

def updateRow(row, toRemove):
    newRow = row
    newRow = update_basic_fields(newRow)
    newRow = hgvsUpdate(newRow)
    newRow["Allele_Frequency"] = selectAlleleFrequency(newRow)
    newRow["Genomic_Coordinate_hg37"] = EMPTY
    newRow["Genomic_Coordinate_hg36"] = EMPTY
    newRow["Source_URL"] = EMPTY
    newRow["Discordant"] = EMPTY
    for item in toRemove:
        del newRow[item]
    return(newRow)

def update_basic_fields(row):
    row["Gene_Symbol"] = row["Gene_symbol(ENIGMA)"]
    row["Genomic_Coordinate"] = row["Genomic_Coordinate(ENIGMA)"]
    row["Reference_Sequence"] = row["Reference_sequence(ENIGMA)"]
    row["HGVS_cDNA"] = row["HGVS_cDNA(ENIGMA)"]
    row["BIC_Identifier"] = row["BIC_Nomenclature(ENIGMA)"]
    row["HGVS_Protein"] = row["HGVS_protein(ENIGMA)"]
    row["Protein_Change"] = row["Abbrev_AA_change(ENIGMA)"]
    return row

def unpackHgvs(string):
    hgvsString = re.sub("\[", "", string)
    hgvsString = re.sub("\]", "", hgvsString)
    hgvsString = re.sub("\'", "", hgvsString)
    firstHgvsString = hgvsString.split("|")[0]
    transcript = firstHgvsString.split(":")[0]
    cleanedHgvsString = re.sub(transcript+":", "", hgvsString)
    return(transcript, cleanedHgvsString)


def hgvsUpdate(row):
    if row["HGVS_cDNA"] == EMPTY:
        if  row["HGVS(ClinVar)"] != EMPTY:
            (row["Reference_Sequence"], 
             row["HGVS_cDNA"]) = unpackHgvs(row["HGVS(ClinVar)"])
        elif  row["HGVS_cDNA(LOVD)"] != EMPTY:
            (row["Reference_Sequence"], 
             row["HGVS_cDNA"]) = unpackHgvs(row["HGVS_cDNA(LOVD)"])
        elif  row["HGVS_cDNA(exLOVD)"] != EMPTY:
            (row["Reference_Sequence"], 
             row["HGVS_cDNA"]) = unpackHgvs(row["HGVS_cDNA(exLOVD)"])
    if row["HGVS_Protein"] == EMPTY:
        if  row["HGVS_protein(LOVD)"] != EMPTY:
            row["HGVS_Protein"] = row["HGVS_protein(LOVD)"]
        elif  row["HGVS_protein(exLOVD)"] != EMPTY:
            row["HGVS_Protein"] = row["HGVS_protein(exLOVD)"]
    return row

def BICUpdate(row):
    if row["BIC_Identifier"] == EMPTY:
        if row["BIC_Nomenclature(BIC)"] != EMPTY:
            return row["BIC_Nomenclature(BIC)"]
        elif row["BIC_Nomenclature(exLOVD)"] != EMPTY:
            return row["BIC_Nomenclature(exLOVD)"]
            

def selectAlleleFrequency(row):
    if row["Allele_frequency(ExAC)"] != EMPTY:
        return row["Allele_frequency(ExAC)"]
    elif row["Minor_allele_frequency(ESP)"] != EMPTY:
        return row["Minor_allele_frequency(ESP)"].split(',')[-1]
    elif row["Allele_frequency(1000_Genomes)"] != EMPTY:
        return row["Allele_frequency(1000_Genomes)"]
    else:
        return EMPTY

if __name__ == "__main__":
    main()

