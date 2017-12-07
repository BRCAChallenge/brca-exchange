#!/usr/bin/env python

'''
calcVarPriors

Parses a tsv file (default built.tsv) containing variant information and for each variant in file 
calculates either the prior probability of pathogenicity or a prior ENGIMA classification based on variant type and variant location
'''

import argparse
import csv
import requests
import sys
import time
import json
import re
from calcMaxEntScanMeanStd import fetch_gene_coordinates

# Here are the canonical BRCA transcripts in ENSEMBL nomenclature
BRCA1_CANONICAL = "ENST00000357654"
BRCA2_CANONICAL = "ENST00000380152"

# Rest Ensembl server
SERVER = "http://rest.ensembl.org"

def checkSequence(sequence):
    '''Checks if a given sequence contains acceptable nucleotides returns True if sequence is comprised entirely of acceptable bases'''
    acceptableBases = ["A", "C", "T", "G", "N", "R", "Y"]
    if len(sequence) > 0:
        for base in sequence:
            if base not in acceptableBases:
                return False
        return True
    else:
        return False


def getVarStrand(variant):
    '''Given a variant, returns the coding strand based on the variant's gene_symbol'''
    varGene = variant["Gene_Symbol"]

    if varGene == "BRCA1": 
        return '-'
    elif varGene == "BRCA2":
        return '+'
    else:
        return ""


def _make_request(url):
    '''Makes request to API and returns json file'''
    req = requests.get(url, headers = {"Content-Type": "application/json"})
    
    if req.status_code == 429 and 'Retry-After' in req.headers:
        retry = float(req.headers['Retry-After'])
        time.sleep(retry)
        return _make_request(url, varData)

    if not req.ok:
        req.raise_for_status()
        sys.exit()

    return req.json()


def getVarConsequences(variant):
    '''
    Given a variant, uses Ensembl VEP API to get variant consequences
    (e.g. intron variant, frameshift variant, missense variant)
    using variant chromosome, Hg38 start, Hg38 end, and alternate allele as input for API
    returns a string detailing consequences of variant
    '''
    ext = "/vep/human/region/"

    # varStrand always 1 because all alternate alleles and positions refer to the plus strand
    varStrand = 1
    varAlt = variant["Alt"]
    
    if variant["Chr"] not in ["13", "17"]:
        return "unable_to_determine"
    else:
        for base in varAlt:
            # API only works for alt alleles that are composed of the 4 canonical bases
            if base not in ["A", "C", "G", "T"]:
                return "unable_to_determine"     
           
        query = "%s:%s-%s:%s/%s?" % (variant["Chr"], variant["Hg38_Start"],
                                     variant["Hg38_End"], varStrand, varAlt)
    
        req_url = SERVER+ext+query
        jsonOutput = _make_request(req_url)
    
        assert(len(jsonOutput) == 1)
        assert(jsonOutput[0].has_key("transcript_consequences"))
        # below is to extract variant consequence from json file
        for gene in jsonOutput[0]["transcript_consequences"]:
            if gene.has_key("transcript_id"):
                # need to filter for canonical BRCA1 transcript
                if re.search(BRCA1_CANONICAL, gene["transcript_id"]):
                    return gene["consequence_terms"][0]
                # need to filter for canonical BRCA2 transcript
                elif re.search(BRCA2_CANONICAL, gene["transcript_id"]):
                    return gene["consequence_terms"][0]
    

def getVarType(variant):
    '''
    Returns a string describing type of variant 
    -substitution, deletion, insertion, delins, other
    depending on variant reference and alternate alleles
    '''
    varRef = variant["Ref"]
    varAlt = variant["Alt"]
    acceptableRefSeq = checkSequence(varRef)
    acceptableAltSeq = checkSequence(varAlt)
    
    if acceptableRefSeq == True and acceptableAltSeq == True: 
        if len(varRef) == len(varAlt):
            if len(varRef) == 1:
                return "substitution"
            else:
                return "delins"
        else:
            # variant is an indel or other variant type
            if len(varRef) > len(varAlt):
                if len(varAlt) == 1:
                    return "deletion"
                else:
                    return "delins"
            elif len(varRef) < len(varAlt):
                if len(varRef) == 1:
                    return "insertion"
                else:
                    return "delins"
            else:
                # variant is not an indel or substitution variant
                return "other"
    else:
        # not acceptable ref seq and alt seq, variant will not be handled by code
        return "other"


def getExonBoundaries(variant):
    '''
    Given a variant, returns the exon boundaries for the variant's transcript in a dictionary with format:
    key = exon number, value = dictionary with exon start and exon end for specific exon
    Uses function implemented in calcMaxEntScanMeanStd to get data for variant's transcript
    '''
    varTranscript = variant["Reference_Sequence"]
    transcriptData = fetch_gene_coordinates(str(varTranscript))

    # parse exon starts and exon ends
    transcriptData["exonStarts"] = re.sub(",(\s)*$", "", transcriptData["exonStarts"])
    transcriptData["exonEnds"] = re.sub(",(\s)*$", "", transcriptData["exonEnds"])
    if transcriptData["strand"] == "+":
        exonStarts = transcriptData["exonStarts"].split(",")
        exonEnds = transcriptData["exonEnds"].split(",")
    else:
        # transcript is on '-' strand
        exonStarts = list(reversed(transcriptData["exonEnds"].split(",")))
        exonEnds = list(reversed(transcriptData["exonStarts"].split(",")))

    varExons = {}
    varExonCount = transcriptData["exonCount"]
    exonCount = 1
    while exonCount <= varExonCount:
        exonStart = int(exonStarts[exonCount - 1])
        exonEnd = int(exonEnds[exonCount - 1])
        if varTranscript == "NM_007294.3":
            if exonCount >= 4:
                # because transcript NM_007294.3 does not include an exon 4, goes from exon 3 to exon 5
                exonName = "exon" + str(exonCount + 1)
            else:
                exonName = "exon" + str(exonCount)
        else:
            exonName = "exon" + str(exonCount)
        varExons[exonName] = {"exonStart": exonStart,
                              "exonEnd": exonEnd}
        exonCount += 1

    return varExons

def getRefSpliceDonorBoundaries(variant):
    varExons = getExonBoundaries(variant)
    varStrand = getVarStrand(variant)
    donorBoundaries = {}
    for exon in varExons.keys():
        exonEnd = int(varExons[exon]["exonEnd"])
        if varStrand == "+":
            donorStart = exonEnd - 3 + 1
            donorEnd = exonEnd + 6
        else:
            # varStrand == "-"
            donorStart = exonEnd + 3
            donorEnd = exonEnd - 6 + 1
        donorBoundaries[exon] = {"donorStart": donorStart,
                                 "donorEnd": donorEnd}

    return donorBoundaries


def getVarLocation(variant):
    '''Given a variant, returns location of variant using Ensembl API for variant annotation'''
    # TO DO - Implement this function using Ensembl API so that variant location is identified
    varLoc = '-'
    return varLoc


def getVarDict(variant):
    '''
    Given input data, returns a dictionary containing information for each variant in input
    Dictionary key is variant HGVS_cDNA and value is a dictionary containing variant gene, variant chromosome, 
    variant strand, variant genomic coordinate, variant type, and variant location
    '''
    varStrand = getVarStrand(variant)
    varType = getVarType(variant)
    varLoc = getVarLocation(variant)

    varDict = {"varGene": variant["Gene_Symbol"],
               "varChrom": variant["Chr"],
               "varStrand": varStrand,
               "varGenCoordinate": variant["Pos"],
               "varType": varType,
               "varLoc": varLoc,
               "varHGVScDNA": variant["pyhgvs_cDNA"]}

    return varDict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputFile", default="built.tsv", help="File with variant information")
    parser.add_argument('-o', "--outputFile", help="File where results will be output")
    parser.add_argument('-b', "--boundaries", default="ENIGMA",
                        help="Specifies which boundaries (ENIGMA or PRIORS) to use for clinically important domains")
    args = parser.parse_args()    
    
    inputData = csv.DictReader(open(args.inputFile, "r"), delimiter="\t")
    for variant in inputData:
        varDict = getVarDict(variant)
    # TO DO - create conditional to account for user selected boundaries
    newColumns = ["varType", "varLoc", "pathProb", "ENIGMAClass", "donorVarMES",
                  "donorVarZ", "donorRefMES", "donorRefZ", "accVarMES", "accVarZ",
                  "accRefMES", "accRefZ", "deNovoMES", "deNovoZ", "spliceSite",
                  "spliceRescue", "frameshift", "CNV", "spliceFlag"]
    # TO DO - create built_with_priors (copy of built) and append new columns
    
if __name__ == "__main__":
    main()
