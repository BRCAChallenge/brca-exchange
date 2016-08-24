#!/usr/bin/env python

import argparse
import csv
import re

EMPTY = "-"
FIELDS_TO_REMOVE=["Gene_symbol_ENIGMA", "Genomic_Coordinate",
                  "Reference_sequence_ENIGMA", "HGVS_cDNA_ENIGMA",
                  "BIC_Nomenclature_ENIGMA", "Abbrev_AA_change_ENIGMA",
                  "HGVS_protein_ENIGMA", "Protein_ClinVar",
                  "HGVS_ClinVar", "HGVS_cDNA_LOVD", "HGVS_cDNA_exLOVD",
                  "HGVS_protein_LOVD", "dna_change_genomic_LOVD",
                  "HGVS_protein_exLOVD",
                  "polyPhen2_result_ESP", 
                  "BIC_Designation_BIC", "BIC_Nomenclature_exLOVD"]
FIELDS_TO_ADD=["Gene_Symbol", "Reference_Sequence",
               "HGVS_cDNA", "Hg38_Start", "Hg38_End",
               "HGVS_RNA", "BIC_Identifier", 
               "HGVS_Protein", 
               "Protein_Change", "Allele_Frequency", 
               "Max_Allele_Frequency",
               "Genomic_Coordinate_hg38",
               "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36", 
               "Source_URL", "Discordant", "Synonyms",
               "Pathogenicity_default", "Pathogenicity_research"]
FIELDS_TO_RENAME={ "Genomic_Coordinate" : "Genomic_Coordinate_Hg38",
                   "Gene_Symbol_ENIGMA" : "Gene_Symbol",
                   "BIC_Nomenclature_ENIGMA" : "BIC_Identifier",
                   "Abbrev_AA_ENIGMA" : "Protein_Change"
                   }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        default="/hive/groups/cgl/brca/release1.0/merged_withVEP_cleaned.csv")
    parser.add_argument("-o", "--output",
                        default="/hive/groups/cgl/brca/release1.0/aggregated.csv")
    args = parser.parse_args()

    csvIn = csv.DictReader(open(args.input, "r"), delimiter='\t')
    outputColumns = setOutputColumns(csvIn.fieldnames, FIELDS_TO_REMOVE,
                                     FIELDS_TO_ADD)
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter='\t',
                            fieldnames=outputColumns)
    csvOut.writerow(dict((fn,fn) for fn in outputColumns))
    rowCount = 0
    for row in csvIn:
        rowCount += 1
        print "working on row", rowCount
        csvOut.writerow(updateRow(row, FIELDS_TO_REMOVE))

        
def setOutputColumns(fields, toRemove, toAdd, toRename):
    newFields = []
    for oldName, newName in toRename:
        newFields.remove(oldName)
        newFields.append(newName)
    for item in fields:
        newFields.append(item)
    for item in toRemove:
        print "removing", item
        newFields.remove(item)
    for item in toAdd:
        newFields.append(item)
    return(newFields)

def updateRow(row, toRemove):
    #if row.has_key(None):
    #    for item in row.keys():
    #        print item, row[item]
    newRow = update_basic_fields(newRow)
    newRow = hgvsUpdate(newRow)
    newRow["Allele_Frequency"] = selectAlleleFrequency(newRow)
    newRow["Max_Allele_Frequency"] = selectMaxAlleleFrequency(newRow)
    newRow["Discordant"] = checkDiscordantStatus(newRow)
    newRow["Source_URL"] = setSourceUrls(newRow)
    newRow["Synonyms"] = setSynonym(row)
    newRow["Genomic_Coordinate_hg37"] = EMPTY
    newRow["Genomic_Coordinate_hg36"] = EMPTY
    for item in toRemove:
        del newRow[item]
    return(newRow)

def update_basic_fields(row):
    row["Genomic_Coordinate_hg38"] = row["Genomic_Coordinate"]
    row["Hg38_Start"] = row['Pos']
    row["Hg38_End"] = int(row["Hg38_Start"]) + len(row["Ref"]) - 1
    row["Gene_Symbol"] = row["Gene_symbol_ENIGMA"]
    if row["Gene_Symbol"] == EMPTY:
        if re.search("^chr17", row["Genomic_Coordinate_hg38"] ):
            row["Gene_Symbol"] = "BRCA1"
        else:
            row["Gene_Symbol"] = "BRCA2"
    row["Reference_Sequence"] = row["Reference_sequence_ENIGMA"]
    row["HGVS_cDNA"] = row["HGVS_cDNA_ENIGMA"]
    row["BIC_Identifier"] = row["BIC_Nomenclature_ENIGMA"]
    row["HGVS_RNA"] = EMPTY
    row["HGVS_Protein"] = row["HGVS_protein_ENIGMA"]
    if row["HGVS_Protein"] == EMPTY:
        row["HGVS_Protein"] = row["Protein_ClinVar"]
    if row["HGVS_Protein"] == EMPTY:
        row["HGVS_Protein"] = row["HGVS_protein_LOVD"]
    if row["HGVS_Protein"] == EMPTY:
        row["HGVS_Protein"] = row["HGVS_protein_exLOVD"]
    row["HGVS_Protein"] = re.sub(row["Reference_Sequence"] + ":", "",
                                 row["HGVS_Protein"])
    row["Protein_Change"] = row["Abbrev_AA_change_ENIGMA"]
    row["Pathogenicity_default"] = row["Clinical_significance_ENIGMA"]
    if row["Pathogenicity_default"] == EMPTY:
        row["Pathogenicity_default"] = "Not Yet Classified"
    if row["Pathogenicity_default"] == "Benign":
        row["Pathogenicity_default"] = "Benign / Little Clinical Significance"
    patho_research = ""
    delimiter = ""
    if row["Clinical_significance_ENIGMA"] != EMPTY:
        patho_research = row["Clinical_significance_ENIGMA"] + "(ENIGMA)"
        delimiter = "; "
    if row["Clinical_Significance_ClinVar"] != EMPTY:
        patho_research = "%s%s%s (ClinVar)" % (patho_research, delimiter,
                                               row["Clinical_Significance_ClinVar"])
        delimiter = "; "
    if row["Clinical_classification_BIC"] != EMPTY:
        patho_research = "%s%s%s (BIC)" % (patho_research, delimiter,
                                           row["Clinical_classification_BIC"])
    row["Pathogenicity_research"] = patho_research
    return row

def unpackHgvs(hgvsString):
    firstHgvsString = hgvsString.split(",")[0]
    if re.search(":", firstHgvsString):
        transcript = firstHgvsString.split(":")[0]
        suffix = re.sub(transcript+":", "", hgvsString)
    elif re.search(".c.", firstHgvsString):
        transcript = firstHgvsString.split(".c.")[0]
        suffix = re.sub(transcript+".", "", hgvsString)
    return(transcript, suffix)


def hgvsUpdate(row):
    if row["HGVS_cDNA"] == EMPTY:
        if  row["HGVS_ClinVar"] != EMPTY:
            (row["Reference_Sequence"], 
             row["HGVS_cDNA"]) = unpackHgvs(row["HGVS_ClinVar"])
        elif  row["HGVS_cDNA_LOVD"] != EMPTY:
            (row["Reference_Sequence"], 
             row["HGVS_cDNA"]) = unpackHgvs(row["HGVS_cDNA_LOVD"])
        elif  row["HGVS_cDNA_exLOVD"] != EMPTY:
            (row["Reference_Sequence"], 
             row["HGVS_cDNA"]) = unpackHgvs(row["HGVS_cDNA_exLOVD"])
    if row["HGVS_Protein"] == EMPTY:
        if  row["HGVS_protein_LOVD"] != EMPTY:
            row["HGVS_Protein"] = row["HGVS_protein_LOVD"]
        elif  row["HGVS_protein_exLOVD"] != EMPTY:
            row["HGVS_Protein"] = row["HGVS_protein_exLOVD"]
    return row

def BICUpdate(row):
    if row["BIC_Identifier"] == EMPTY:
        if row["BIC_Nomenclature_BIC"] != EMPTY:
            return row["BIC_Designation_BIC"]
        elif row["BIC_Nomenclature_exLOVD"] != EMPTY:
            return row["BIC_Nomenclature_exLOVD"]
            

def selectAlleleFrequency(row):
    if row["Allele_frequency_ExAC"] != EMPTY:
        return "%s (ExAC" % row["Allele_frequency_ExAC"]
    elif row["Minor_allele_frequency_ESP"] != EMPTY:
        return "%s (ESP)" % row["Minor_allele_frequency_ESP"].split(',')[-1]
    elif row["Allele_frequency_1000_Genomes"] != EMPTY:
        return "%s (1000 Genomes)" % row["Allele_frequency_1000_Genomes"]
    else:
        return EMPTY

def selectMaxAlleleFrequency(newRow):
    maxFreq = 0
    maxFreqString = EMPTY
    if newRow["Minor_allele_frequency_ESP"] != EMPTY:
        #print newRow["Minor_allele_frequency_ESP"]
        tokens = newRow["Minor_allele_frequency_ESP"].split(",")
        if len(tokens) >= 1:
            ea = float(tokens[0]) / 100
            if ea > maxFreq:
                maxFreq = ea
                maxFreqString = "%f (EA from ESP)" % ea
        if len(tokens) > 2:
            aa = float(tokens[1]) / 100
            if aa >= maxFreq:
                maxFreq = aa
                maxFreqString = "%f (AA from ESP)" % aa
    if newRow["EUR_Allele_frequency_1000_Genomes"] != EMPTY:
        eur_af = float(newRow["EUR_Allele_frequency_1000_Genomes"])
        if eur_af > maxFreq:
            maxFreq = eur_af
            maxFreqString = "%f (EUR from 1000 Genomes)" % eur_af
    if newRow["AFR_Allele_frequency_1000_Genomes"] != EMPTY:
        afr_af = float(newRow["AFR_Allele_frequency_1000_Genomes"])
        if afr_af > maxFreq:
            maxFreq = afr_af
            maxFreqString = "%f (AFR from 1000 Genomes)" % afr_af
    if newRow["AMR_Allele_frequency_1000_Genomes"] != EMPTY:
        amr_af = float(newRow["AMR_Allele_frequency_1000_Genomes"])
        if amr_af > maxFreq:
            maxFreq = amr_af
            maxFreqString = "%f (AMR from 1000 Genomes)" % amr_af
    if newRow["EAS_Allele_frequency_1000_Genomes"] != EMPTY:
        eas_af = float(newRow["EAS_Allele_frequency_1000_Genomes"])
        if eas_af > maxFreq:
            maxFreq = eas_af
            maxFreqString = "%f (EAS from 1000 Genomes)" % eas_af
    if newRow["SAS_Allele_frequency_1000_Genomes"] != EMPTY:
        sas_af = float(newRow["SAS_Allele_frequency_1000_Genomes"])
        if sas_af > maxFreq:
            maxFreq = sas_af
            maxFreqString = "%f (SAS from 1000 Genomes)" % sas_af
    return(maxFreqString)


def checkDiscordantStatus(row):
    hasPathogenicClassification = False
    hasBenignClassification = False
    for column in (row["Clinical_Significance_ClinVar"], row["Clinical_significance_ENIGMA"]):
        for item in column.split(","):
            if re.search("^pathogenic$", item.lower()):
                hasPathogenicClassification = True
            if re.search("^pathologic$", item.lower()):
                hasPathogenicClassification = True
            if re.search("^likely_pathogenic$", item.lower()):
                hasPathogenicClassification = True
            if re.search("^probable_pathogenic$", item.lower()):
                hasPathogenicClassification = True
            if re.search("^benign$", item.lower()):
                hasBenignClassification = True
            if re.search("^probably_not_pathogenic$", item.lower()):
                hasBenignClassification = True
            if re.search("^likely_benign$", item.lower()):
                hasBenignClassification = True
            if re.search("^no_known_pathogenicity$", item.lower()):
                hasBenignClassification = True
            if re.search("^variant_of_unknown_significance$", item.lower()):
                hasBenignClassification = True
            if re.search("^uncertain_significance$", item.lower()):
                hasBenignClassification = True
    for item in row["Clinical_classification_BIC"].split(","):
        if re.search("class 5", item.lower()):
            hasPathogenicClassification = True
        if re.search("class 1", item.lower()):
            hasBenignClassification = True
    if hasPathogenicClassification and hasBenignClassification:
        return "Discordant"
    else:
        return "Concordant"

def setSourceUrls(row):
    url = ""
    delimiter = ""
    if row["URL_ENIGMA"] != EMPTY:
        for thisURL in row["URL_ENIGMA"].split(','):
            url = "%s%s%s" % (url, delimiter, thisURL)
            delimiter = ", "
    if row["SCV_ClinVar"] != EMPTY:
        for thisSCV in row["SCV_ClinVar"].split(','):
            variantUrl = "http://www.ncbi.nlm.nih.gov/clinvar/?term="+ thisSCV
            url = "%s%s%s" % (url, delimiter, variantUrl)
            delimiter=", "
    if url != "":
        return url
    else:
        return EMPTY

def setSynonym(row):
    delimiter = ""
    synonym = ""
    if row["BIC_Nomenclature_ENIGMA"] != EMPTY:
        for thisBic in row["BIC_Nomenclature_ENIGMA"].split(','):
            synonym = "%s%s%s" % (synonym, delimiter, thisBic)
            delimiter = ","
    if row["BIC_Nomenclature_exLOVD"] != EMPTY:
        for thisBic in row["BIC_Nomenclature_exLOVD"].split(','):
            synonym = "%s%s%s" % (synonym, delimiter, thisBic)
            delimiter = ","
    if row["BIC_Designation_BIC"] != EMPTY:
        for thisBic in row["BIC_Designation_BIC"].split(','):
            synonym = "%s%s%s" % (synonym, delimiter, thisBic)
            delimiter = ","
    return synonym


if __name__ == "__main__":
    main()

