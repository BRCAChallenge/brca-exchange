#!/usr/bin/env python

import argparse
import csv
import re

EMPTY = "-"
FIELDS_TO_REMOVE=["Gene_symbol(ENIGMA)", "Genomic_Coordinate",
                  "Reference_sequence(ENIGMA)", "HGVS_cDNA(ENIGMA)",
                  "BIC_Nomenclature(ENIGMA)", "Abbrev_AA_change(ENIGMA)",
                  "HGVS_protein(ENIGMA)", 
                  "HGVS(ClinVar)", "HGVS_cDNA(LOVD)", "HGVS_cDNA(exLOVD)",
                  "HGVS_protein(LOVD)", "HGVS_protein(exLOVD)",
                  "polyPhen2_result(ESP)", 
                  "BIC_Designation(BIC)", "BIC_Nomenclature(exLOVD)"]
FIELDS_TO_ADD=["Gene_Symbol", "Reference_Sequence",
               "HGVS_cDNA", "BIC_Identifier", "HGVS_Protein", 
               "Protein_Change", "Allele_Frequency", 
               "Genomic_Coordinate_hg38",
               "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36", 
               "Source_URL", "Discordant", "Other_HGVS_cDNA",
               "Pathogenicity_default", "Pathogenicity_research"]

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
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter='\t',
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

def cleanColumns(row):
    """For each column, remove any pickling that might be delimiting the
    column, to generate proper text"""
    for column in row.keys():
        value = row[column]
        value = re.sub("^\[", "", value)
        value = re.sub("\]$", "", value)
        value = re.sub("^'", "", value)
        value = re.sub("'$", "", value)
        value = re.sub("'\]\|\['", ",", value)
        value = re.sub("\]\|\[", ",", value)
        value = re.sub("\|", ",", value)
        row[column] = value
    return row

def updateRow(row, toRemove):
    newRow = cleanColumns(row)
    newRow = update_basic_fields(newRow)
    newRow = hgvsUpdate(newRow)
    newRow["Allele_Frequency"] = selectAlleleFrequency(newRow)
    newRow["Discordant"] = checkDiscordantStatus(newRow)
    newRow["Genomic_Coordinate_hg37"] = EMPTY
    newRow["Genomic_Coordinate_hg36"] = EMPTY
    newRow["Source_URL"] = EMPTY
    for item in toRemove:
        del newRow[item]
    return(newRow)

def update_basic_fields(row):
    row["Genomic_Coordinate_hg38"] = row["Genomic_Coordinate"]
    row["Gene_Symbol"] = row["Gene_symbol(ENIGMA)"]
    if row["Gene_Symbol"] == EMPTY:
        if re.search("^chr17", row["Genomic_Coordinate_hg38"] ):
            row["Gene_Symbol"] = "BRCA1"
        else:
            row["Gene_Symbol"] = "BRCA2"
    row["Reference_Sequence"] = row["Reference_sequence(ENIGMA)"]
    row["Genomic_Coordinate_hg38"] = row["Genomic_Coordinate"]
    row["HGVS_cDNA"] = row["HGVS_cDNA(ENIGMA)"]
    row["BIC_Identifier"] = row["BIC_Nomenclature(ENIGMA)"]
    row["HGVS_Protein"] = row["HGVS_protein(ENIGMA)"]
    if row["HGVS_Protein"] == EMPTY:
        row["HGVS_Protein"] = row["HGVS_protein(LOVD)"]
    if row["HGVS_Protein"] == EMPTY:
        row["HGVS_Protein"] = row["HGVS_protein(exLOVD)"]
    row["HGVS_Protein"] = re.sub(row["Reference_Sequence"] + ":", "",
                                 row["HGVS_Protein"])
    row["Protein_Change"] = row["Abbrev_AA_change(ENIGMA)"]
    row["Pathogenicity_default"] = row["Clinical_significance(ENIGMA)"]
    if row["Pathogenicity_default"] == EMPTY:
        row["Pathogenicity_default"] = "Not Yet Classified"
    patho_research = ""
    delimiter = ""
    if row["Clinical_significance(ENIGMA)"] != EMPTY:
        patho_research = row["Clinical_significance(ENIGMA)"] + "(ENIGMA)"
        delimiter = "; "
    if row["Clinical_Significance(ClinVar)"] != EMPTY:
        patho_research = "%s%s%s (ClinVar)" % (patho_research, delimiter,
                                               row["Clinical_Significance(ClinVar)"])
        delimiter = "; "
    if row["Clinical_classification(BIC)"] != EMPTY:
        patho_research = "%s%s%s (BIC)" % (patho_research, delimiter,
                                           row["Clinical_classification(BIC)"])
    row["Pathogenicity_research"] = patho_research
    if len(row["SIFT(VEP)"]) == 0:
        row["SIFT(VEP)"] = EMPTY
    if len(row["PolyPhen(VEP)"]) == 0:
        row["PolyPhen(VEP)"] = EMPTY
    return row

def unpackHgvs(hgvsString):
    firstHgvsString = hgvsString.split("|")[0]
    if re.search(":", firstHgvsString):
        transcript = firstHgvsString.split(":")[0]
        suffix = re.sub(transcript+":", "", hgvsString)
    elif re.search(".c.", firstHgvsString):
        transcript = firstHgvsString.split(".c.")[0]
        suffix = re.sub(transcript+".", "", hgvsString)
    return(transcript, suffix)


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
            return row["BIC_Designation(BIC)"]
        elif row["BIC_Nomenclature(exLOVD)"] != EMPTY:
            return row["BIC_Nomenclature(exLOVD)"]
            

def selectAlleleFrequency(row):
    if row["Allele_frequency(ExAC)"] != EMPTY:
        return "%s (ExAC)" % row["Allele_frequency(ExAC)"]
    elif row["Minor_allele_frequency(ESP)"] != EMPTY:
        return "%s (ESP)" % row["Minor_allele_frequency(ESP)"].split(',')[-1]
    elif row["Allele_frequency(1000_Genomes)"] != EMPTY:
        return "%s (1000 Genomes)" % row["Allele_frequency(1000_Genomes)"]
    else:
        return EMPTY


def checkDiscordantStatus(row):
    hasPathogenicClassification = False
    if re.search("pathogenic", row["Clinical_Significance(ClinVar)"].lower()):
        hasPathogenicClassification = True
    if re.search("pathogenic", row["Clinical_significance(ENIGMA)"].lower()):
        hasPathogenicClassification = True
    if re.search("pathologic", row["Clinical_Significance(ClinVar)"].lower()):
        hasPathogenicClassification = True
    if re.search("class 5", row["Clinical_classification(BIC)"].lower()):
        hasPathogenicClassification = True
    hasBenignClassification = False
    if re.search("benign", row["Clinical_Significance(ClinVar)"].lower()):
        hasBenignClassification = True
    if re.search("no_known_pathogenicity", 
                 row["Clinical_Significance(ClinVar)"].lower()):
        hasBenignClassification = True
    if re.search("benign", row["Clinical_significance(ENIGMA)"].lower()):
        hasBenignClassification = True
    if re.search("class 1", row["Clinical_classification(BIC)"].lower()):
        hasBenignClassification = True
    if hasPathogenicClassification and hasBenignClassification:
        return "True"
    else:
        return "False"


if __name__ == "__main__":
    main()

