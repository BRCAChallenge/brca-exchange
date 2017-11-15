#!/usr/bin/env python

'''
calcVarPriors

Parses a tsv file (default built.tsv) containing variant information and for each variant in file 
calculates either the prior probability of pathogenicity or a prior ENGIMA classification based on variant type and variant location
'''

import argparse
import csv
import re
from calcMaxEntScanMeanStd import fetch_gene_coordinates


def checkSequence(sequence):
    '''Checks if a given sequence contains acceptable nucleotides returns True if sequence is comprised entirely of acceptable bases'''
    acceptableBases = ["A", "C", "T", "G", "N", "R", "Y"]
    badBases = 0
    if len(sequence) > 0:
        for base in sequence:
            if base not in acceptableBases:
                badBases += 1
            if badBases == 0:
                acceptableSequence = True
            else:
                # badBases > 0
                acceptableSequence = False
    else:
        # len(sequence) = 0
        acceptableSequence = False
    return acceptableSequence


def getVarType(variant):
    '''
    Returns a string describing type of variant (substitution, deletion, insertion, delins, other) depending on variant reference and alternate alleles
    '''
    acceptableRefSeq = checkSequence(variant["Ref"])
    acceptableAltSeq = checkSequence(variant["Alt"])
    if acceptableRefSeq == True and acceptableAltSeq == True: 
        if len(variant["Ref"]) == len(variant["Alt"]):
            if len(variant["Ref"]) == 1:
                varType = "substitution"
            else:
                varType = "delins"
        else:
            # variant is an indel or other variant type
            if len(variant["Ref"]) > len(variant["Alt"]):
                if len(variant["Alt"]) == 1:
                    varType = "deletion"
                else:
                    # delins type variant
                    varType = "delins"
            elif len(variant["Ref"]) < len(variant["Alt"]):
                if len(variant["Ref"]) == 1:
                    varType = "insertion"
                else:
                    # delins type variant
                    varType = "delins"
            else:
                # variant is not an indel or substitution variant
                varType = "other"
    else:
        # not acceptable ref seq and alt seq, variant will not be handled by code
        varType = "other"
    return varType


def getExonBoundaries(transcript):
    '''
    Given a transcript, returns a dictionary containing the exon number as the key and a tuple containing the exon start and exon end as the key for each exon
    Uses fetch_gene_coordinates from calcMaxEntScanMeanStd to get information for a particular transcript
    '''
    varTranscript = fetch_gene_coordinates(str(transcript))
    # parse exon starts and exon ends
    varTranscript['exonStarts'] = re.sub(",(\s)*$", "", varTranscript['exonStarts'])
    varTranscript['exonEnds'] = re.sub(",(\s)*$", "", varTranscript['exonEnds'])
    if varTranscript['strand'] == "+":
        exonStarts = varTranscript['exonStarts'].split(',')
        exonEnds = varTranscript['exonEnds'].split(',')
    else:
        # variant is on '-' strand
        exonStarts = list(reversed(varTranscript['exonEnds'].split(',')))
        exonEnds = list(reversed(varTranscript['exonStarts'].split(',')))
    varExons = {}
    varExonCount = varTranscript['exonCount']
    exonCount = 1
    while exonCount <= varExonCount:
        exonStart = int(exonStarts[exonCount-1])
        exonEnd = int(exonEnds[exonCount-1])
        exonName = "exon" + str(exonCount)
        varExons[exonName] = {"exonStart":exonStart,
                              "exonEnd":exonEnd}
        exonCount += 1
    return varExons

def varOutsideBoundaries(variant):
    '''
    Given a variant, determines if variant genomic position is outside transcript boundaries, handles variants on positive and negative strand
    Returns either True (outside boundaries) or False (inside boundaries)
    '''
    varGenPos = int(variant['Pos'])
    varRefTranscript = variant["Reference_Sequence"]
    varTranscript = fetch_gene_coordinates(str(varRefTranscript))
    varStrand = varTranscript['strand']
    outsideBoundaries = False
    if varStrand == '+':
        txnStart = int(varTranscript['txStart'])
        txnEnd = int(varTranscript['txEnd'])
        if varGenPos < txnStart or varGenPos > txnEnd:
            outsideBoundaries = True
    else:
        # varStrand == '-'
        txnStart = int(varTranscript['txEnd'])
        txnEnd = int(varTranscript['txStart'])
        if varGenPos > txnStart or varGenPos < txnEnd:
            outsideBoundaries = True
    return outsideBoundaries
    
def varInExon(variant):
    '''
    Given a variant, determines if variant genomic position is in an exon, handles variants on positive and negative strand
    Returns either True (in an exon) or False (not in an exon)
    '''
    varGenPos = int(variant["Pos"])
    varRefTranscript = variant["Reference_Sequence"]
    varExons = getExonBoundaries(str(varRefTranscript))
    varTranscript = fetch_gene_coordinates(str(varRefTranscript))
    varStrand = varTranscript["strand"]
    inExon = False
    for exon in varExons:
        exonBounds = varExons[exon]
        exonStart = int(exonBounds['exonStart'])
        exonEnd = int(exonBounds['exonEnd'])
        if varStrand == '+':
            if varGenPos >= exonStart and varGenPos <= exonEnd:
                inExon = True
        else:
            # varStrand == '-'
            if varGenPos <= exonStart and varGenPos >= exonEnd:
                inExon = True
    return inExon


def getVarLocation(variant):
    '''Given a variant, returns general location of variant (outside transcript boundaries, in an exon, or an intron)'''
    varOutsideBounds = varOutsideBoundaries(variant)
    varExon = varInExon(variant)
    if varOutsideBounds == True:
        varLoc = "outsideBoundaries"
    else:
        # varOutsideBoundaries == False
        if varExon == True:
            varLoc = "inExon"
        else:
            # inExon == False
            varLoc = "inIntron"
    return varLoc


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
    varLoc = getVarLocation(variant)
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
    
    newColumns = ["varType", "varLoc", "pathProb", "ENIGMAClass", "donorVarMES", "donorVarZ", "donorRefMES", "donorRefZ", "accVarMES", "accVarZ", "accRefMES", "accRefZ", "deNovoMES", "deNovoZ", "spliceSite", "spliceRescue", "frameshift", "CNV", "spliceFlag"]
    # TO DO - create built_with_priors (copy of built) and append new columns
    
if __name__ == "__main__":
    main()
