#!/usr/bin/env python
"""add VEP result to merged data"""
import argparse
import csv
import json
import pandas as pd
import re
import requests
import sys

# Here are the canonical BRCA transcripts in ENSEMBL nomenclature
BRCA1_CANONICAL = "ENST00000357654"
BRCA2_CANONICAL = "ENST00000380152"

VEP_TRANSCRIPT_CONSEQUENCES = {
    "Sift_Score" : "sift_score",
    "Sift_Prediction" : "sift_prediction",
    "Polyphen_Score" : "polyphen_score", 
    "Polyphen_Prediction" : "polyphen_prediction"
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        default="/hive/groups/cgl/brca/release1.0/merged.tsv")
    parser.add_argument("-o", "--output",
                        default="/hive/groups/cgl/brca/release1.0/merged_withVEP_cleaned.tsv")
    parser.add_argument("-v", "--vep",  help="VEP output, run in advance for all variants",
                        default="/cluster/home/mollyzhang/release1.0/data/VEP/vep_output_3_3_2016.vcf")
    args = parser.parse_args()
    csvIn = csv.DictReader(open(args.input, "r"), delimiter='\t')
    outputColumns = setOutputColumns(csvIn.fieldnames, VEP_TRANSCRIPT_CONSEQUENCES)
    csvOut = csv.DictWriter(open(args.output, "w"), delimiter='\t',
                            fieldnames=outputColumns)
    csvOut.writerow(dict((fn,fn) for fn in outputColumns))
    rowCount = 0
    for row in csvIn:
        rowCount += 1
        row = addVepResults(row, VEP_TRANSCRIPT_CONSEQUENCES)
        csvOut.writerow(row)

def setOutputColumns(fields, toAdd):
    newFields = []
    for item in fields:
        newFields.append(item)
    for item in toAdd:
        newFields.append(item)
    return(newFields)

def addVepResults(row, vepTranscriptConsequenceFields):
    # Initialize to the default output
    defaultOutput = "-"
    for label, field in vepTranscriptConsequenceFields.iteritems():
        row[label] = defaultOutput
    # Assemble a query.  As of this writing, the API doesn't seem to work for 
    # anything other than simple missense substitutions.  Skip this process
    # unless both the reference and alt alleles are one of the four canonical 
    # bases, and are different from each other.
    if row["Ref"] in ['A', 'C', 'G', 'T'] \
            and row['Alt'] in ['A', 'C', 'G', 'T'] \
            and row["Ref"] != row["Alt"]:
        server = "http://rest.ensembl.org"
        ext = "/vep/human/hgvs/"
        hgvs = "%s:g.%s:%s>%s?" % (row["Chr"], row["Pos"], 
                                   row["Ref"], row["Alt"])
        req = requests.get(server+ext+hgvs, 
                           headers={ "Content-Type" : "application/json"})
        if not req.ok:
            req.raise_for_status()
            sys.exit()
        jsonOutput = req.json()
        assert(len(jsonOutput) == 1)
        assert(jsonOutput[0].has_key("transcript_consequences"))
        correctEntry = None
        for entryThisGene in jsonOutput[0]["transcript_consequences"]:
            if entryThisGene.has_key("transcript_id"):
                if re.search(BRCA1_CANONICAL, entryThisGene["transcript_id"]):
                    correctEntry = entryThisGene
                elif re.search(BRCA2_CANONICAL, 
                               entryThisGene["transcript_id"]):
                    correctEntry = entryThisGene
        for label, field in vepTranscriptConsequenceFields.iteritems():
            if correctEntry != None:
                if correctEntry.has_key(field):
                    row[label] = correctEntry[field]
    return(row)

if __name__ == "__main__":
    main()
