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
    revCompDict = {"A": "T",
                   "T": "A",
                   "C": "G",
                   "G": "C",
                   "N": "N"}

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
    
    varParsed =  {"typeHGVS": varHGVS.kind,
                  "varRef": varHGVS.ref_allele,
                  "varAlt": varHGVS.alt_allele,
                  "varType": varHGVS.mutation_type}
    
    return varParsed


def getVarInfo(variant):
    '''
    Given a variant, returns a dictionary containing: 
    HGVS_cDNA, varType, varGenRef, varGenAlt, and cDNA ref and alt alleles
    If varGene is BRCA1, cDNA ref and alt are the reverse complement 
    If varGene is BRCA2, cDNA ref and alt are normal
    '''
    varGene = variant["Gene_Symbol"]
    varHGVS = variant["pyhgvs_cDNA"]
    varParse = parseVar(varHGVS)

    varInfo = {"varGenRef": variant["Ref"],
               "varGenAlt": variant["Alt"],
               "varcDNAHGVS": varHGVS,
               "varType": varParse["varType"],
               "varGene": varGene}

    if varGene == "BRCA2":
        cDNARef = varParse["varRef"]
        cDNAAlt = varParse["varAlt"]
    elif varGene == "BRCA1":
        # gets reverse complement because BRCA1 on minus strand
        cDNARef = getRevComp(varParse["varRef"])
        cDNAAlt = getRevComp(varParse["varAlt"])
    else:
        # varGene is not BRCA1 or BRCA2
        cDNARef = varParse["varRef"]
        cDNAAlt = varParse["varAlt"]
        
    varInfo["cDNARef"] = cDNARef
    varInfo["cDNAAlt"] = cDNAAlt

    return varInfo


def checkVarHGVS(variant):
    '''
    Checks that ref and alt alleles for a variant are consistent with parsed HGVS_cDNA
    If consistent, returns True
    Otherwise returns a list with information about the variant that includes:
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
            varEqual = False
    elif varType == "del":
        # compares cDNA ref allele with second base to end of the genomic ref allele
        # first base of genomic ref allele is unchaged between genomic ref and alt alleles
        # and HGVS_cDNA does not include unchanged bases
        # so is excluded from comparison between cDNA ref and genomic ref alleles
        tempVarGenRef = varGenRef[1:]
        varEqual = tempVarGenRef == cDNARef
    elif varType == "dup":
        # compares cDNA ref allele with second base to end of the genomic alt allele
        # first base of genomic alt allele is unchanged between genomic ref and alt alleles
        # so it is excluded from comparison
        # compare cDNA ref to genomic alt because cDNA alt is actually duplicated cDNA ref sequence
        tempVarGenAlt = varGenAlt[1:]
        varEqual = tempVarGenAlt == cDNARef
    else:
        # compares cDNA alt allele with second base to end of genomic alt allele
        # first base of genomic alt allele is unchanged between genomic ref and alt alleles
        # and HGVS_cDNA does not include unchanged bases
        # so first based is excluded from comparison between cDNA alt and genomic alt alleles
        tempVarGenAlt = varGenAlt[1:]
        varEqual = tempVarGenAlt == cDNAAlt
        
    if varEqual == True:
        return varEqual
    else:
        varData = [varInfo["varGene"],
                   varInfo["varcDNAHGVS"],
                   varType,
                   varGenRef,
                   varGenAlt,
                   cDNARef,
                   cDNAAlt]
        
        return varData
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputFile", help = "File with variant HGVS to be checked")
    parser.add_argument('-o', "--outputFile", help = "Output file that contains variants with inconsistent genomic and cDNA HGVS")
    args = parser.parse_args()

    inputData = csv.DictReader(open(args.inputFile, "r"), delimiter="\t")
    outputData = csv.writer(open(args.outputFile, "w"), delimiter="\t")

    headers = ["Gene_Symbol", "HGVS_cDNA", "varType", "genRef", "genAlt", "cDNARef", "cDNAAlt"]
    outputData.writerow(headers)

    inconVars = 0
    for variant in inputData:
        checkVar = checkVarHGVS(variant)
        if checkVar == True:
            pass
        else:
            inconVars += 1
            outputData.writerow(checkVar)
    
    print "Number of inconsitent genomic and cDNA HGVS variants is:" + str(inconVars)

    outputData.close()
            
if __name__ == "__main__":
    main()
