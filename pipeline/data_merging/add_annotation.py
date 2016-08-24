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
    "SIFT_SCORE" : "sift_score",
    "SIFT_PREDICTION" : "sift_prediction",
    "POLYPHEN_SCORE" : "polyphen_score", 
    "POLYPHEN_PREDICTION" : "polyphen_prediction"
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
    # anything other than simple missense substitutions.  Skip this process unless
    # both the reference and alt alleles are one of the four canonical bases.
    if row["Ref"] in ['A', 'C', 'G', 'T'] and row['Alt'] in ['A', 'C', 'G', 'T'] \
            and row["Ref"] != row["Alt"]:
        server = "http://rest.ensembl.org"
        ext = "/vep/human/hgvs/"
        hgvs = "%s:g.%s:%s>%s?" % (row["Chr"], row["Pos"], row["Ref"], row["Alt"])
        print "workign on", hgvs
        req = requests.get(server+ext+hgvs, headers={ "Content-Type" : "application/json"})
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
                elif re.search(BRCA2_CANONICAL, entryThisGene["transcript_id"]):
                    correctEntry = entryThisGene
        for label, field in vepTranscriptConsequenceFields.iteritems():
            if correctEntry != None:
                if correctEntry.has_key(field):
                    row[label] = correctEntry[field]
    return(row)

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
        items = line.strip().split("\t")
        if line_num == 1:
            vep_fields = [i+"_VEP" for i in VEP_FIELDS]
            items += vep_fields
        else:
            genome_coor = items[2]
            additional_items = []
            coordinate_this_variant = ""
            for coordinate in genome_coor.split(","):
                if coordinate in vep_result_dict.keys():
                    coordinate_this_variant = coordinate
            if coordinate_this_variant == "":
                for column in VEP_FIELDS:
                    additional_items.append("-")
            else:
                for column in VEP_FIELDS:
                    this_cell = vep_result_dict[coordinate_this_variant][column]
                    additional_items.append(this_cell)
            items += additional_items
        new_line = "\t".join(items) + "\n"
        f_out.write(new_line)



def save_VEP_to_dict(vepFile):
    print("processing VEP output")
    f = open(vepFile, "r")
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
    #print "hello world"
    main()
