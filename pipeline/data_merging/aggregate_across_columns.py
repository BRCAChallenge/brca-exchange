#!/usr/bin/env python
import argparse
import copy
import csv
import re
import utilities

EMPTY = "-"
FIELDS_TO_REMOVE = ["Protein_ClinVar",
                    "Description_ClinVar",
                    "Summary_Evidence_ClinVar",
                    "Review_Status_ClinVar",
                    "HGVS_ClinVar",
                    "HGVS_cDNA_exLOVD",
                    "HGVS_protein_LOVD",
                    "HGVS_protein_exLOVD",
                    "polyPhen2_result_ESP",
                    "BIC_Designation_BIC",
                    "BIC_Nomenclature_exLOVD"]
FIELDS_TO_ADD = ["Hg38_Start", "Hg38_End", "Hg37_Start", "Hg37_End",
                 "Hg36_Start", "Hg36_End",
                 "HGVS_RNA",
                 "Allele_Frequency",
                 "Max_Allele_Frequency",
                 "Genomic_Coordinate_hg37", "Genomic_Coordinate_hg36",
                 "Source_URL", "Discordant", "Synonyms",
                 "Pathogenicity_expert", "Pathogenicity_all"]
FIELDS_TO_RENAME = {"Gene_symbol_ENIGMA": "Gene_Symbol",
                    "Genomic_Coordinate": "Genomic_Coordinate_hg38",
                    "Reference_sequence_ENIGMA": "Reference_Sequence",
                    "Abbrev_AA_change_ENIGMA": "Protein_Change",
                    "HGVS_cDNA_ENIGMA": "HGVS_cDNA",
                    "HGVS_protein_ENIGMA": "HGVS_Protein",
                    "BIC_Nomenclature_ENIGMA": "BIC_Nomenclature"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        default="/hive/groups/cgl/brca/release1.0/merged_withVEP_cleaned.csv")
    parser.add_argument("-o", "--output",
                        default="/hive/groups/cgl/brca/release1.0/aggregated.csv")
    args = parser.parse_args()

    csvIn = csv.DictReader(open(args.input, "r"), delimiter='\t')
    outputColumns = setOutputColumns(csvIn.fieldnames, FIELDS_TO_REMOVE,
                                     FIELDS_TO_ADD, FIELDS_TO_RENAME)
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter='\t',
                            fieldnames=outputColumns)
    csvOut.writerow(dict((fn, fn) for fn in outputColumns))
    rowCount = 0
    for row in csvIn:
        rowCount += 1
        csvOut.writerow(updateRow(row, FIELDS_TO_RENAME, FIELDS_TO_REMOVE))
    print "Process complete, aggregated %s variants." % (rowCount)


def setOutputColumns(fields, toRemove, toAdd, toRename):
    newFields = []
    for item in fields:
        newFields.append(item)
    for oldName, newName in toRename.iteritems():
        newFields.remove(oldName)
        newFields.append(newName)
    for item in toRemove:
        newFields.remove(item)
    for item in toAdd:
        newFields.append(item)
    return(newFields)


def updateRow(row, toRename, toRemove):
    newRow = copy.deepcopy(row)
    newRow = update_basic_fields(row, toRename)
    (newRow["Reference_Sequence"], newRow["HGVS_cDNA"]) = hgvsCdnaUpdate(newRow)
    newRow["HGVS_Protein"] = hgvsProteinUpdate(row)
    newRow["BIC_Nomenclature"] = BICUpdate(row)
    (newRow["Pathogenicity_expert"],
     newRow["Pathogenicity_all"]) = pathogenicityUpdate(newRow)
    newRow["Allele_Frequency"] = selectAlleleFrequency(newRow)
    newRow["Max_Allele_Frequency"] = selectMaxAlleleFrequency(newRow)
    newRow["Discordant"] = checkDiscordantStatus(newRow)
    newRow["Source_URL"] = setSourceUrls(newRow)
    newRow["Synonyms"] = setSynonym(row)
    newRow["Genomic_Coordinate_hg37"] = EMPTY
    newRow["Genomic_Coordinate_hg36"] = EMPTY
    newRow["Hg37_Start"] = EMPTY
    newRow["Hg37_End"] = EMPTY
    newRow["Hg36_Start"] = EMPTY
    newRow["Hg36_End"] = EMPTY
    for item in toRemove:
        del newRow[item]
    for item in newRow:
        assert item != None
        assert len(item) > 0
    return(newRow)


def update_basic_fields(row, columnsToReplace):
    for key, value in columnsToReplace.iteritems():
        row[value] = row[key]
        del row[key]
    row["Hg38_Start"] = row['Pos']
    row["Hg38_End"] = int(row["Hg38_Start"]) + len(row["Ref"]) - 1
    if row["Gene_Symbol"] == EMPTY:
        if re.search("^chr17", row["Genomic_Coordinate_hg38"]):
            row["Gene_Symbol"] = "BRCA1"
        else:
            row["Gene_Symbol"] = "BRCA2"
    row["HGVS_RNA"] = EMPTY
    return row


def unpackHgvs(hgvsString):
    firstHgvsString = hgvsString.split(",")[0]
    if re.search(":", firstHgvsString):
        transcript = firstHgvsString.split(":")[0]
        suffix = re.sub(transcript+":", "", hgvsString)
    elif re.search(".c.", firstHgvsString):
        transcript = firstHgvsString.split(".c.")[0]
        suffix = re.sub(transcript+".", "", hgvsString)
    elif re.search(".g.", firstHgvsString):
        transcript = firstHgvsString.split(".g.")[0]
        suffix = re.sub(transcript+".", "", hgvsString)
    #TODO: Could use an else case, transcript assignment errors could occurs here.
    return(transcript, suffix)


def hgvsCdnaUpdate(row):
    refSequence = row["Reference_Sequence"]
    hgvs = row["HGVS_cDNA"]
    if hgvs == EMPTY:
        if row["HGVS_ClinVar"] != EMPTY:
            (refSequence, hgvs) = unpackHgvs(row["HGVS_ClinVar"])
        elif row["HGVS_cDNA_LOVD"] != EMPTY:
            (refSequence, hgvs) = unpackHgvs(row["HGVS_cDNA_LOVD"])
        elif row["HGVS_cDNA_exLOVD"] != EMPTY:
            (refSequence, hgvs) = unpackHgvs(row["HGVS_cDNA_exLOVD"])
    return(refSequence, hgvs)


def hgvsProteinUpdate(row):
    protein = row["HGVS_Protein"]
    if protein == EMPTY:
        if row["Protein_ClinVar"] == EMPTY:
            protein = row["Protein_ClinVar"]
        elif row["HGVS_protein_LOVD"] != EMPTY:
            protein = row["HGVS_protein_LOVD"]
        elif row["HGVS_protein_exLOVD"] != EMPTY:
            protein = row["HGVS_protein_exLOVD"]
    # 8/24/16: this is an error condition that should not occur.  There should be
    # an assertion checking the input that no value is empty.  We know in practice
    # that there are cases where this assertion fails, currently.  I'm fixing the
    # input code.  When that fix is in, this "if" can go away: protein will
    # always be a string with length of at least 1.
    if protein != None:
        protein = re.sub(row["Reference_Sequence"] + ":", "",
                         protein)
    return protein


def BICUpdate(row):
    bic = row["BIC_Nomenclature"]
    bic = re.sub("\|", ",", bic)
    if bic == EMPTY:
        if row["BIC_Designation_BIC"] != EMPTY:
            bic = row["BIC_Designation_BIC"]
        elif row["BIC_Nomenclature_exLOVD"] != EMPTY:
            bic = row["BIC_Nomenclature_exLOVD"]
    return bic


def pathogenicityUpdate(row):
    pathoExpert = row["Clinical_significance_ENIGMA"]
    if pathoExpert == EMPTY:
        pathoExpert = "Not Yet Reviewed"
    if pathoExpert == "Benign":
        pathoExpert = "Benign / Little Clinical Significance"
    pathoAll = ""
    delimiter = ""
    if row["Clinical_significance_ENIGMA"] != EMPTY:
        pathoAll = row["Clinical_significance_ENIGMA"] + "(ENIGMA)"
        delimiter = "; "
    if row["Clinical_Significance_ClinVar"] != EMPTY:
        pathoAll = "%s%s%s (ClinVar)" % (pathoAll, delimiter,
                                               row["Clinical_Significance_ClinVar"])
        delimiter = "; "
    if row["Clinical_classification_BIC"] != EMPTY:
        pathoAll = "%s%s%s (BIC)" % (pathoAll, delimiter,
                                           row["Clinical_classification_BIC"])
    if pathoAll == "":
        pathAll = EMPTY
    return(pathoExpert, pathoAll)


def selectAlleleFrequency(row):
    if row["Allele_frequency_ExAC"] != EMPTY:
        return "%s (ExAC minus TCGA)" % row["Allele_frequency_ExAC"]
    elif row["Minor_allele_frequency_percent_ESP"] != EMPTY:
        # Percent must be converted to a fraction
        return "%s (ESP)" % (float(row["Minor_allele_frequency_percent_ESP"].split(',')[-1])/100)
    elif row["Allele_frequency_1000_Genomes"] != EMPTY:
        return "%s (1000 Genomes)" % row["Allele_frequency_1000_Genomes"]
    else:
        return EMPTY


def selectMaxAlleleFrequency(newRow):
    # MaxAF was removed from the UI and is now automatically set to '-' as of 8/11/17.
    return '-'


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
            variantUrl = "http://www.ncbi.nlm.nih.gov/clinvar/?term=" + thisSCV
            url = "%s%s%s" % (url, delimiter, variantUrl)
            delimiter = ", "
    if url != "":
        return url
    else:
        return EMPTY


def setSynonym(row):
    delimiter = ""
    synonym = ""
    if row["BIC_Nomenclature"] != EMPTY:
        for thisBic in row["BIC_Nomenclature"].split(','):
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
    return (synonym)


if __name__ == "__main__":
    main()
