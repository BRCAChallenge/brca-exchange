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

# clinically important domain boundaries
brca1CIDomains = {"enigma": {"ring": {"domStart": 43124096,
                                      "domEnd": 43104260},
                             "brct": {"domStart": 43070966,
                                      "domEnd": 43045705}},
                  "priors": {"initiation": {"domStart": 43124096,
                                            "domEnd": 43124094},
                             "ring": {"domStart": 43124084,
                                      "domEnd": 43104875},
                             "brct": {"domStart": 43070966,
                                      "domEnd": 43045705}}}

brca2CIDomains = {"enigma": {"dnb": {"domStart": 32356433,
                                     "domEnd": 32396954}},
                  "priors": {"initiation": {"domStart": 32316461,
                                            "domEnd": 32316463},
                             "palb2": {"domStart": 32316491,
                                       "domEnd": 32319108},
                             "dnb": {"domStart": 32356433,
                                     "domEnd": 32396954},
                             "tr2/rad5": {"domStart": 32398318,
                                          "domEnd": 32398428}}}

# BRCA1/BRCA2 grey zone boundaries
greyZones = {"BRCA2": {"greyZoneStart": 32398438,
                       "greyZoneEnd": 32398488}}


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

def checkWithinBoundaries(varStrand, varGenPos, boundaryStart, boundaryEnd):
    '''
    Checks whether a position (varGenPos) is within certain boundaries (boundaryStart and boundaryEnd)
    Dependent on the variant's transcript strand (varStrand)
    Function is inclusive of boundaryStart and boundaryEnd
    '''
    if varStrand == "+":
        if varGenPos >= boundaryStart and varGenPos <= boundaryEnd:
            return True
    elif varStrand == "-":
        if varGenPos <= boundaryStart and varGenPos >= boundaryEnd:
            return True
    else:
        return False

def varOutsideBoundaries(variant):
    '''Given a variant, determines if variant is outside transcript boundaries'''
    varGenPos = int(variant["Pos"])
    varTranscript = variant["Reference_Sequence"]
    transcriptData = fetch_gene_coordinates(str(varTranscript))
    varStrand = getVarStrand(variant)
    if varStrand == "+":
        txnStart = int(transcriptData["txStart"])
        txnEnd = int(transcriptData["txEnd"])
        if varGenPos < txnStart or varGenPos > txnEnd:
            return True
    else:
        txnStart = int(transcriptData["txEnd"])
        txnEnd = int(transcriptData["txStart"])
        if varGenPos > txnStart or varGenPos < txnEnd:
            return True

    return False

def varInUTR(variant):
    '''
    Given a variant, if variant is inside transcript boundaries, 
    determines if variant is in 3' or 5' UTR of transcript
    '''
    varOutBounds = varOutsideBoundaries(variant)
    if varOutBounds == False:
        varGenPos = int(variant["Pos"])
        varTranscript = variant["Reference_Sequence"]
        transcriptData = fetch_gene_coordinates(varTranscript)
        varStrand = getVarStrand(variant)
        if varStrand == "+":
            tsnStart = int(transcriptData["cdsStart"])
            tsnEnd = int(transcriptData["cdsEnd"])
            if varGenPos < tsnStart:
                return True
            elif varGenPos > tsnEnd:
                return True
        else:
            tsnStart = int(transcriptData["cdsEnd"])
            tsnEnd = int(transcriptData["cdsStart"])
            if varGenPos > tsnStart:
                return True
            elif varGenPos < tsnEnd:
                return True
    return False
                
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
    # exonCount starts at 1 because there is no exon 0 in transcripts
    # and exonCount is directly used for exon name
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
    '''
    Given a variant, returns the splice donor boundaries 
    (splice donor region is last 3 bases in exon and frist 6 bases in intron) 
    for the variant's transcript in a dictionary with the format:
    key = exon number, value = dictionary with donor start and donor end for exon
    '''
    varExons = getExonBoundaries(variant)
    if variant["Gene_Symbol"] == "BRCA1":
        del varExons["exon24"]
    elif variant["Gene_Symbol"] == "BRCA2":
        del varExons["exon27"]
    varStrand = getVarStrand(variant)
    donorBoundaries = {}
    for exon in varExons.keys():
        exonEnd = int(varExons[exon]["exonEnd"])
        if varStrand == "+":
            # - 3 + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, donor start is 5' to exon end for + strand transcripts
            donorStart = exonEnd - 3 + 1
            donorEnd = exonEnd + 6
        else:
            donorStart = exonEnd + 3
            # - 6 + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, donor end is 5' to exon end for - strand transcripts
            donorEnd = exonEnd - 6 + 1
        donorBoundaries[exon] = {"donorStart": donorStart,
                                 "donorEnd": donorEnd}

    return donorBoundaries

def getRefSpliceAcceptorBoundaries(variant):
    '''
    Given a variant, returns the splice acceptor boundaries
    (splice acceptor region is 20 bases before exon and first 3 bases in exon)
    for the variant's transcript in a dictionary with the format:
    key = exon number, value = a dictionary with acceptor start and acceptor end for exon
    '''
    varExons = getExonBoundaries(variant)
    if variant["Gene_Symbol"] == "BRCA1" or variant["Gene_Symbol"] == "BRCA2":
        del varExons["exon1"]
    varStrand = getVarStrand(variant)
    acceptorBoundaries = {}
    for exon in varExons.keys():
        exonStart = int(varExons[exon]["exonStart"])
        if varStrand == "+":
            # -20 + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, acceptor start is 5' to exon start for + strand transcripts
            acceptorStart = exonStart - 20 + 1
            acceptorEnd = exonStart + 3
        else:            
            acceptorStart = exonStart + 20
            # -3 + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, acceptor end is 5' to exon start for - strand transcripts
            acceptorEnd = exonStart - 3 + 1
        acceptorBoundaries[exon] = {"acceptorStart": acceptorStart,
                                    "acceptorEnd": acceptorEnd}

    return acceptorBoundaries

def varInExon(variant):
    '''
    Given a variant, determines if variant genomic position is inside transcript boundaries
    AND if variant is in an exon
    Returns true if variant is in an exon
    '''
    varOutBounds = varOutsideBoundaries(variant)
    if varOutBounds == False:
        varGenPos = int(variant["Pos"])
        varExons = getExonBoundaries(variant)
        varStrand = getVarStrand(variant)
        for exon in varExons.keys():
            exonStart = varExons[exon]["exonStart"]
            exonEnd = varExons[exon]["exonEnd"]
            if varStrand == "+":
                if varGenPos > exonStart and varGenPos <= exonEnd:
                    return True
            else:
                withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, exonStart, exonEnd)
                if withinBoundaries == True:
                    return True
    return False
    
    
def varInSpliceDonor(variant):
    '''
    Given a variant, determines if a variant is in reference transcript's splice donor region
    splice donor region = last 3 bases in exon and first 6 bases in intron
    Returns True if variant in splice donor region
    '''
    varGenPos = int(variant["Pos"])
    varStrand = getVarStrand(variant)
    donorBounds = getRefSpliceDonorBoundaries(variant)
    for exon in donorBounds.keys():
        donorStart = donorBounds[exon]["donorStart"]
        donorEnd = donorBounds[exon]["donorEnd"]
        withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, donorStart, donorEnd)
        if withinBoundaries == True:
            return True
    return False

def varInSpliceAcceptor(variant):
    '''
    Given a variant, determines if a variant is in reference transcript's splice acceptor region
    splice acceptor region = 20 bases preceding exon and first 3 bases in exon
    Returns True if variant in splice acceptor region
    '''
    varGenPos = int(variant["Pos"])
    varStrand = getVarStrand(variant)
    acceptorBounds = getRefSpliceAcceptorBoundaries(variant)
    for exon in acceptorBounds.keys():
        acceptorStart = acceptorBounds[exon]["acceptorStart"]
        acceptorEnd = acceptorBounds[exon]["acceptorEnd"]
        withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, acceptorStart, acceptorEnd)
        if withinBoundaries == True:
            return True
    return False

def varInCiDomain(variant, boundaries):
    '''
    Given a variant, determines if variant is in a clinically important domain
    Second argument determiens which boundaries (ENIGMA or PRIORS) are used for CI domains
    Returns True if variant in CI domain
    '''
    varGenPos = int(variant["Pos"])
    varGene = variant["Gene_Symbol"]
    varStrand = getVarStrand(variant)
    inExon = varInExon(variant)
    if inExon == True:
        if varGene == "BRCA1":
            for domain in brca1CIDomains[boundaries].keys():
                domainStart = brca1CIDomains[boundaries][domain]["domStart"]
                domainEnd = brca1CIDomains[boundaries][domain]["domEnd"]
                withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, domainStart, domainEnd)
                if withinBoundaries == True:
                    return True
        elif varGene == "BRCA2":
            for domain in brca2CIDomains[boundaries].keys():
                domainStart = brca2CIDomains[boundaries][domain]["domStart"]
                domainEnd = brca2CIDomains[boundaries][domain]["domEnd"]
                withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, domainStart, domainEnd)
                if withinBoundaries == True:
                    return True
    return False
                
def varInGreyZone(variant):
    '''
    Given a variant, determines if variant is in the grey zone
    Returns True if variant is in the grey zone
    Returns False if variant is NOT in the grey zone
    '''
    varGenPos = int(variant["Pos"])
    varGene = variant["Gene_Symbol"]
    varStrand = getVarStrand(variant)
    if varGene == "BRCA2":
        greyZoneStart = greyZones[varGene]["greyZoneStart"]
        greyZoneEnd = greyZones[varGene]["greyZoneEnd"]
        withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, greyZoneStart, greyZoneEnd)
        if withinBoundaries == True:
            return True
    return False

def varAfterGreyZone(variant):
    '''
    Given a variant, determines if variant is after the gene grey zone
    Returns True if variant is after the grey zone
    '''
    varGenPos = int(variant["Pos"])
    varGene = variant["Gene_Symbol"]
    # checks that varGene == BRCA2 because only BRCA2 has a grey zone
    if varGene == "BRCA2":
        inUTR = varInUTR(variant)
        inGreyZone = varInGreyZone(variant)
        # makes sure that variant is not in grey zone
        # and makes sure that variant is not in UTR of gene
        if inGreyZone == False and inUTR == False:
            greyZoneEnd = greyZones[varGene]["greyZoneEnd"]
            if varGenPos > greyZoneEnd:
                return True
    return False
                
def getVarLocation(variant, boundaries):
    '''
    Given a variant, returns the variant location as below
    Second argument is for CI domain boundaries (PRIORS or ENIGMA)
    '''
    varOutBounds = varOutsideBoundaries(variant)
    if varOutBounds == True:
        return "outside_transcript_boundaries_variant"
    inExon = varInExon(variant)
    inSpliceDonor = varInSpliceDonor(variant)
    inSpliceAcceptor = varInSpliceAcceptor(variant)
    if inExon == True:
        inCiDomain = varInCiDomain(variant, boundaries)
        if inCiDomain == True and inSpliceDonor == True:
            return "CI_splice_donor_variant"
        if inCiDomain == True and inSpliceAcceptor == True:
            return "CI_splice_acceptor_variant"
        if inCiDomain == True:
            return "CI_domain_variant"
        if inSpliceDonor == True:
            return "splice_donor_variant"
        if inSpliceAcceptor == True:
            return "splice_acceptor_variant"
        inGreyZone = varInGreyZone(variant)
        if inGreyZone == True:
            return "grey_zone_variant"
        afterGreyZone = varAfterGreyZone(variant)
        if afterGreyZone == True:
            return "after_grey_zone_variant"
        inUTR = varInUTR(variant)
        if inUTR == True:
            return "UTR_variant"
        return "exon_variant"
    else:    
        if inSpliceDonor == True:
            return "splice_donor_variant"
        if inSpliceAcceptor == True:
            return "splice_acceptor_variant"
        inUTR = varInUTR(variant)
        if inUTR == True:
            return "UTR_variant"
        return "intron_variant"
    
def getVarDict(variant, boundaries):
    '''
    Given input data, returns a dictionary containing information for each variant in input
    Dictionary key is variant HGVS_cDNA and value is a dictionary containing variant gene, variant chromosome, 
    variant strand, variant genomic coordinate, variant type, and variant location
    '''
    varStrand = getVarStrand(variant)
    varType = getVarType(variant)
    varLoc = getVarLocation(variant, boundaries)

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
