#!/usr/bin/env python

import argparse
import csv
from Bio.Seq import Seq
import pyhgvs as hgvs

'''
checkVarHGVS

Takes a tsv as input (usually built.tsv) and checks if the variant ref and alt alleles are consistent with the parsed HGVS_cDNA string for that variant

Prints the total number of variants with inconsistencies

Outputs a file that contains information for each inconsistent variant
'''

def getRevComp(sequence):
    '''Given a sequence returns the reverse complement'''
    revCompDict = {"A":"T",
                   "T":"A",
                   "C":"G",
                   "G":"C",
                   "N":"N"}
    if len(sequence) == 1:
        revComp = revCompDict[sequence]
    else:
        revComp = str(Seq(sequence).reverse_complement())
    return revComp

def parseVar(variantHGVS):
    '''
    Parses the given variant HGVS and returns a dictionary containing: 
    HGVS type, variant type, ref allele, and alt allele
    ''' 
    varHGVS = hgvs.HGVSName(str(variantHGVS))
    typeHGVS = varHGVS.kind
    varRef = varHGVS.ref_allele
    varAlt = varHGVS.alt_allele
    varType = varHGVS.mutation_type
    varParsed =  {"typeHGVS":typeHGVS,
                  "varRef":varRef,
                  "varAlt":varAlt,
                  "varType":varType}
    return varParsed

def getVarInfo(variant):
    '''
    Given a variant, returns a dictionary containing HGVS_cDNA, varType, varGenRef, varGenAlt, and depending on if varGene is BRCA1 or BRCA2 either the reverse complement or normal cDNA ref and alt alleles respectively
    ''' 
    varGenRef = variant["Ref"]
    varGenAlt = variant["Alt"]
    varcDNAHGVS = variant["pyhgvs_cDNA"]
    varGene = variant["Gene_Symbol"]
    varParse = parseVar(varcDNAHGVS)
    varType = varParse["varType"]
    varInfo = {"varGenRef":varGenRef,
               "varGenAlt":varGenAlt,
               "varcDNAHGVS":varcDNAHGVS,
               "varType":varType,
               "varGene":varGene}
    if varGene == "BRCA2":
        cDNARef = varParse["varRef"]
        cDNAAlt = varParse["varAlt"]
    else:
        # varGene == "BRCA1"
        # gets reverse complement because BRCA1 on minus strand
        cDNARef = getRevComp(varParse["varRef"])
        cDNAAlt = getRevComp(varParse["varAlt"])
    varInfo["cDNARef"] = cDNARef
    varInfo["cDNAAlt"] = cDNAAlt
    return varInfo
    
def checkVarHGVS(variant):
    '''
    Checks that ref and alt alleles for a variant are consistent with parsed HGVS_cDNA
    If consistent, returns True
    Otherwise returns a dictionary with 2 entries
    "Consistent":False
    and
    "varData": a list of information about the variant that includes:
        HGVS_cDNA, variant type, ref allele, alt allele, parsed cDNA ref allele, and parsed cDNA alt allele
    '''
    varInfo = getVarInfo(variant)
    varType = varInfo["varType"]
    varGenRef = varInfo["varGenRef"]
    varGenAlt = varInfo["varGenAlt"]
    cDNARef = varInfo["cDNARef"]
    cDNAAlt = varInfo["cDNAAlt"]
    if varType == ">" or varType == "delins":
        if varGenRef == cDNARef and varGenAlt == cDNAAlt:
            varEqual = True
        else:
            # varGenRef != cDNARef and varGenAlt != cDNAAlt
            varEqual = False
    elif varType == "del":
        tempVarGenRef = varGenRef[1:]
        if tempVarGenRef == cDNARef:
            varEqual = True
        else:
            varEqual = False
    elif varType == "dup":
        tempVarGenAlt = varGenAlt[1:]
        if tempVarGenAlt == cDNARef:
            varEqual = True
        else:
            varEqual = False
    else:
        # varType == "ins":
        tempVarGenAlt = varGenAlt[1:]
        if tempVarGenAlt == cDNAAlt:
            varEqual = True
        else:
            varEqual = False
    if varEqual == True:
        return varEqual
    else:
        # varEqual == False
        varData = varData = [varInfo["varcDNAHGVS"], varType, varGenRef, varGenAlt, cDNARef, cDNAAlt]
        return {"Consistent":varEqual,
                "varData":varData}
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputFile", help = "File with variant HGVS to be checked")
    parser.add_argument('-o', "--outputFile", help = "Output file that contains variants with inconsistent genomic and cDNA HGVS")
    args = parser.parse_args()

    inputData = csv.DictReader(open(args.inputFile, "r"), delimiter="\t")
    outputData = csv.writer(open(args.outputFile, "w"), delimiter="\t")

    headers = ["HGVS_cDNA", "varType", "genRef", "genAlt", "cDNARef", "cDNAAlt"]
    outputData.writerow(headers)

    inconVars = 0
    for variant in inputData:
        checkVar = checkVarHGVS(variant)
        if checkVar == True:
            pass
        else:
            inconVars += 1
            outputData.writerow(checkVar["varData"])

    print "Number of inconsitent genomic and cDNA HGVS variants is:" + str(inconVars)
            
if __name__ == "__main__":
    main()
