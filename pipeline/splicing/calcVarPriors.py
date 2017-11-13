#!/usr/bin/env python

'''
calcVarPriors

Parses a tsv file (default built.tsv) containing variant information and for each variant in file 
calculates either the prior probability of pathogenicity or a prior ENGIMA classification based on variant type and variant location
'''

import argparse
import csv

def getVarType(variant):
    '''
    Returns a string describing type of variant (substitution, deletion, insertion, or delins) depending on variant reference and alternate alleles
    '''
    if len(variant["Ref"]) == len(variant["Alt"]):
        if len(variant["Ref"]) == 1:
            varType = "substitution"
        else:
            varType = "delins"
    else:
        # variant is an indel
        if len(variant["Ref"]) > len(variant["Alt"]):
            if len(variant["Alt"]) == 1:
                varType = "deletion"
            else:
                # delins type variant
                varType = "delins"
        else:
            # len(variant["Ref"]) < len(variant["Alt"])
            if len(variant["Ref"]) == 1:
                varType = "insertion"
            else:
                # delins type variant
                varType = "delins"
    return varType

def getVarDict(variant):
    '''
    Given input data, returns a dictionary containing information for each variant in input
    Dictionary key is variant HGVS_cDNA and value is a dictionary containing variant gene, variant chromosome, 
    variant strand, variant genomic coordinate, variant type, and variant location
    '''
    varHGVS = variant["pyhgvs_cDNA"]
    varGene = variant["Gene_Symbol"]
    varChrom = variant["Chr"]
    if varGene == "BRCA1":
        varStrand = '-'
    else:
        # varGene == "BRCA2"
        varStrand = '+'
    varGenCoordinate = variant["Pos"]
    varType = getVarType(variant)
    varLoc = "-" # until implement varLocation function
    # TO DO - implement varLocation function
    # VarLoc = varLocation(variant)
    varDict = {"varGene":varGene,
               "varChrom":varChrom,
               "varStrand":varStrand,
               "varGenCoordinate":varGenCoordinate,
               "varType":varType,
               "varLoc":varLoc,
               "varHGVScDNA":varHGVS}
    return varDict
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputFile", default="built.tsv", help="File with variant information")
    parser.add_argument('-o', "--outputFile", help="File where results will be output")
    args = parser.parse_args()
    
    inputData = csv.DictReader(open(args.inputFile, "r"), delimiter="\t")
    for variant in inputData:
        varDict = getVarDict(variant)
    
    newColumns = ["varType", "varLoc", "pathProb", "ENIGMAClass", "donorVarMES", "donorVarZ", "donorRefMES", "donorRefZ", "accVarMES", "accVarZ", "accRefMES", "accRefZ", "deNovoMES", "deNovoZ", "spliceSite", "spliceRescue", "frameshift", "CNV"]
    # TO DO - create built_with_priors (copy of built) and append new columns
    
if __name__ == "__main__":
    main()
