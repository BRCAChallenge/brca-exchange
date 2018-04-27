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
import subprocess
import tempfile
import os
from Bio.Seq import Seq
from calcMaxEntScanMeanStd import fetch_gene_coordinates, runMaxEntScan

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

# standard window sizes for splice donor and acceptor
# 3/5/18, defined by Sean Tavtigian and Michael Parsons
# also found in MaxEntScan score splice site definitions (Yeo and Burge 2004)
STD_DONOR_SIZE = 9
STD_DONOR_INTRONIC_LENGTH = 6
STD_DONOR_EXONIC_LENGTH = 3
STD_ACC_SIZE = 23
STD_ACC_INTRONIC_LENGTH = 20
STD_ACC_EXONIC_LENGTH = 3

# standard exonic portion size and de novo acceptor length
# 3/5/18, defined by Sean Tavtigian and Michael Parsons
# stdExonicPortion also found in MaxEntScan score splice site definitions (Yeo and Burge 2004)
STD_EXONIC_PORTION = 3
STD_DE_NOVO_LENGTH = 10

# standard de novo offset
# subject to change if values above change
# default value as of 3/5/18 is 7
STD_DE_NOVO_OFFSET = STD_DE_NOVO_LENGTH - STD_EXONIC_PORTION

# Canonical BRCA transcripts in RefSeq nomenclature
BRCA1_RefSeq = "NM_007294.3"
BRCA2_RefSeq = "NM_000059.3"

# Fetch transcript data for BRCA1/BRCA2 RefSeq transcripts
brca1TranscriptData = fetch_gene_coordinates(BRCA1_RefSeq)
brca2TranscriptData = fetch_gene_coordinates(BRCA2_RefSeq)

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

def getVarChrom(variant):
    '''Given a variant, returns the chromosome based on the variant's gene_symbol'''
    varGene = variant["Gene_Symbol"]

    if varGene == "BRCA1":
        return "chr17"
    elif varGene == "BRCA2":
        return "chr13"
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

def getTranscriptData(referenceSequence):
    '''
    Given a reference sequence (e.g. "NM_007294.3"),
    Returns transcript data for that reference sequencee
    '''
    if referenceSequence == BRCA1_RefSeq:
        return brca1TranscriptData
    elif referenceSequence == BRCA2_RefSeq:
        return brca2TranscriptData
    
def varOutsideBoundaries(variant):
    '''Given a variant, determines if variant is outside transcript boundaries'''
    varGenPos = int(variant["Pos"])
    varTranscript = variant["Reference_Sequence"]
    transcriptData = getTranscriptData(varTranscript)
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
        transcriptData = getTranscriptData(varTranscript)
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
    transcriptData = getTranscriptData(varTranscript)
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

def getRefSpliceDonorBoundaries(variant, intronicLength, exonicLength):
    '''
    Given a variant, intronicLength and exonicLength returns the splice donor boundaries
    intronicLength = number of bp in intron that will be considered as part of splice donor region
    exonicLength = number of bp in exon that will be considered as part of splice donor region
    splice region is the last exonicLength bp in the exon and first intronicLength bp in the intron
    for the variant's transcript in a dictionary with the format:
    key = exon number, value = dictionary with donor start and donor end for exon
    '''
    varExons = getExonBoundaries(variant)
    donorExons = varExons.copy()
    if variant["Gene_Symbol"] == "BRCA1":
        del donorExons["exon24"]
    elif variant["Gene_Symbol"] == "BRCA2":
        del donorExons["exon27"]
    varStrand = getVarStrand(variant)
    donorBoundaries = {}
    for exon in donorExons.keys():
        exonEnd = int(donorExons[exon]["exonEnd"])
        if varStrand == "+":
            # exonicLength + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, donor start is 5' to exon end for + strand transcripts
            donorStart = exonEnd - exonicLength + 1
            donorEnd = exonEnd + intronicLength
        else:
            donorStart = exonEnd + exonicLength
            # intronicLength + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, donor end is 5' to exon end for - strand transcripts
            donorEnd = exonEnd - intronicLength + 1
        donorBoundaries[exon] = {"donorStart": donorStart,
                                 "donorEnd": donorEnd}

    return donorBoundaries

def getSpliceAcceptorBoundaries(variant, intronicLength, exonicLength):
    '''
    Given a variant, intronicLength and exonicLength returns the splice acceptor boundaries
    intronicLength = number of bp in intron that will be considered as part of splice acceptor region
    exonicLength = number of bp in exon that will be considered as part of splice acceptor region
    splice rgion is the last intronicLength bp in the exon and first exonicLength bp in the exon
    for the variant's transcript in a dictionary with the format:
    key = exon number, value = a dictionary with acceptor start and acceptor end for exon
    '''
    varExons = getExonBoundaries(variant)
    acceptorExons = varExons.copy()
    if variant["Gene_Symbol"] == "BRCA1" or variant["Gene_Symbol"] == "BRCA2":
        del acceptorExons["exon1"]
    varStrand = getVarStrand(variant)
    acceptorBoundaries = {}
    for exon in acceptorExons.keys():
        exonStart = int(acceptorExons[exon]["exonStart"])
        if varStrand == "+":
            # intronicLength + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, acceptor start is 5' to exon start for + strand transcripts
            acceptorStart = exonStart - intronicLength + 1
            acceptorEnd = exonStart + exonicLength
        else:
            acceptorStart = exonStart + intronicLength
            # exonicLength + 1 because genomic position in RefSeq starts to the right of the first base
            # which affects 5' side of sequence, acceptor end is 5' to exon start for - strand transcripts
            acceptorEnd = exonStart - exonicLength + 1
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
            exonStart = int(varExons[exon]["exonStart"])
            exonEnd = int(varExons[exon]["exonEnd"])
            if varStrand == "+":
                if varGenPos > exonStart and varGenPos <= exonEnd:
                    return True
            else:
                withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, exonStart, exonEnd)
                if withinBoundaries == True:
                    return True
    return False

def getVarExonNumberSNS(variant):
    '''
    Given a SNS variant, checks that variant is in an exon
    If variant in an exon, returns the number of the exon variant is located within in format "exonN"
    '''
    if varInExon(variant) == True:
        varGenPos = int(variant["Pos"])
        varExons = getExonBoundaries(variant)
        varStrand = getVarStrand(variant)
        for exon in varExons.keys():
            exonStart = varExons[exon]["exonStart"]
            exonEnd = varExons[exon]["exonEnd"]
            if varStrand == "+":
                if varGenPos > exonStart and varGenPos <= exonEnd:
                    return exon
            else:
                withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, exonStart, exonEnd)
                if withinBoundaries == True:
                    return exon

def varInSpliceRegion(variant, donor=False, deNovo=False):
    '''
    Given a variant, determines if a variant is in reference transcript's splice donor/acceptor region
    If donor=True and deNovo=False, checks if variant is in a reference splice donor region
    If donor=True and deNovo=True, checks if variant is in a de novo splice donor region
    If donor=False and deNovo=False, checks if variant is in a reference splice acceptor region
    If donor=False and deNovo=True, checks if variant is in a de novo splice acceptor region
    Returns True if variant is in a splice region, false otherwise
    '''
    if donor == False and deNovo == False:
        regionBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
    elif donor == False and deNovo == True:
        regionBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_DE_NOVO_LENGTH)
    elif donor == True:
        # gets reference donor splice boundaries, if deNovo = True then entireity of exon will be included below
        regionBounds = getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
    for exon in regionBounds.keys():
        if donor == False:
            regionStart = regionBounds[exon]["acceptorStart"]
            regionEnd = regionBounds[exon]["acceptorEnd"]
        else:
            regionStart = regionBounds[exon]["donorStart"]
            regionEnd = regionBounds[exon]["donorEnd"]
        withinBoundaries = checkWithinBoundaries(getVarStrand(variant), int(variant["Pos"]), regionStart, regionEnd)
        if withinBoundaries == True and donor == False:
            return True
        elif donor == True and deNovo == False and withinBoundaries == True:
            return True
        # because de novo donor region includes reference splice donor region and entirity of exon
        elif donor == True and deNovo == True and (withinBoundaries == True or varInExon(variant) == True):
            return True
    return False

def getVarSpliceRegionBounds(variant, donor=False, deNovo=False):
    '''
    Given a variant, checks if variant is in a splice donor/acceptor region
    If donor=True, checks if variant is in a splice donor region and returns boundaries for splice donor region
      *function CANNOT be used to return de novo donor splice region bounds*
    If donor=False and deNovo=False, checks if variant is in a ref splice acceptor region and returns boundaries for splice acceptor region
    If donor=False and deNovo=True, checks if variant is in a de novo splice acceptor region and returns boundaries for that region
    If variant is in a splice region, returns a dictionary with region boundaries where variant is located
    '''
    if varInSpliceRegion(variant, donor=donor, deNovo=deNovo):
        if donor == False:
            if deNovo == False:
                regionBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
            else:
                regionBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_DE_NOVO_LENGTH)
            regionStartKey = "acceptorStart"
            regionEndKey = "acceptorEnd"
        else:        
            regionBounds = getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
            regionStartKey = "donorStart"
            regionEndKey = "donorEnd"
        for exon in regionBounds.keys():
            regionStart = regionBounds[exon][regionStartKey]
            regionEnd = regionBounds[exon][regionEndKey]
            withinBoundaries = checkWithinBoundaries(getVarStrand(variant), int(variant["Pos"]), regionStart, regionEnd)
            if withinBoundaries == True:
                return {"exonName": exon,
                        regionStartKey: regionStart,
                        regionEndKey: regionEnd}    
                
def varInCIDomain(variant, boundaries):
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
    inSpliceDonor = varInSpliceRegion(variant, donor=True, deNovo=False)
    inSpliceAcceptor = varInSpliceRegion(variant, donor=False, deNovo=False)
    if inExon == True:
        inCIDomain = varInCIDomain(variant, boundaries)
        if inCIDomain == True and inSpliceDonor == True:
            return "CI_splice_donor_variant"
        if inCIDomain == True and inSpliceAcceptor == True:
            return "CI_splice_acceptor_variant"
        if inCIDomain == True:
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

def getFastaSeq(chrom, rangeStart, rangeStop, plusSeq=True):
    '''
    Given chromosome (in format 'chr13'), region genomic start position, and
    region genomic end position:
    Returns a string containing the sequence inclusive of rangeStart and rangeStop
    If plusSeq=True, returns plus strand sequence
    If plusSeq=False, returns minus strand sequence
    '''
    if rangeStart < rangeStop:
        regionStart = rangeStart
        regionEnd = rangeStop
    else:
        regionStart = rangeStop
        regionEnd = rangeStart
    
    url = "http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=%s:%d,%d" % (chrom, regionStart, regionEnd)
    req = requests.get(url)
    
    if req.status_code == 429 and 'Retry-After' in req.headers:
        retry = float(req.headers['Retry-After'])
        time.sleep(retry)
        req = requests.get(url)
    
    lines = req.content.split('\n')
    # because sequence is located at index 5 in dictionary
    sequence = lines[5]
    if plusSeq == True:
        return sequence.upper()
    else:
        return str(Seq(sequence.upper()).reverse_complement())

def getSeqLocDict(chrom, strand, rangeStart, rangeStop):
    '''
    Given chromosome, strand, region genomic start position, and region genomic end position
    returns a dictionary containing the genomic position as the key and reference allele as the value
    For minus strand gene (rangeStart > rangeStop), for plus strand gene (rangeStart < rangeStop)
    Always returns plus strand sequence 
    '''
    seqLocDict = {}
    if strand == "-":
        regionStart = int(rangeStop)
        regionEnd = int(rangeStart)
    else:
        regionStart = int(rangeStart)
        regionEnd = int(rangeStop)
    sequence = getFastaSeq(chrom, regionStart, regionEnd, plusSeq=True)
    genPos = regionStart   
    while genPos <= regionEnd:
        for base in sequence:
            seqLocDict[genPos] = base
            genPos += 1
    return seqLocDict

def getAltSeqDict(variant, seqLocDict):
    '''
    Given a variant and a dictionary containing a sequence with bases and their locations,
    returns the dictionary with the alternate allele in place of the reference allele
    at the variant's genomic position
    '''
    varRef = variant["Ref"]
    varAlt = variant["Alt"]
    varGenPos = int(variant["Pos"])
    seqLocDictRef = seqLocDict[varGenPos]
    if varRef == seqLocDictRef:
        altSeqDict = seqLocDict.copy()
        altSeqDict[varGenPos] = varAlt
    return altSeqDict

def getAltSeq(altSeqDict, strand):
    '''
    Given a dictionary containing an alternate sequence with bases and their locations
    and the strand that the alternate allele is on
    Returns a string of the sequence containing the alternate allele
    '''
    sequence = ""
    # to ensure that items in dictionary are sorted numerically
    for key, value in sorted(altSeqDict.items()):
        sequence += altSeqDict[key]
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    return sequence

def getRefAltSeqs(variant, rangeStart, rangeStop):
    '''
    Given a variant, rangeStart, and rangeStop:
    Returns a dicitonary with ref and alt seq for the specified variant and range
    '''
    varChrom = getVarChrom(variant)
    varStrand = getVarStrand(variant)
    if varStrand == "-":
        refSeq = getFastaSeq(varChrom, rangeStart, rangeStop, plusSeq=False)
    else:
        refSeq = getFastaSeq(varChrom, rangeStart, rangeStop, plusSeq=True)
    refSeqDict = getSeqLocDict(varChrom, varStrand, rangeStart, rangeStop)
    altSeqDict = getAltSeqDict(variant, refSeqDict)
    altSeq = getAltSeq(altSeqDict, varStrand)
    return {"refSeq": refSeq,
            "altSeq": altSeq}

def getVarSeqIndexSNS(refSeq, altSeq):
    '''
    Given a reference sequence and alternate sequence for a substitution variant
    Determines the index at which the alt seq differs from the ref seq
      - using zero indexing
    Returns the index at which altSeq differs from refSeq
    For example if the refSeq is ACTG and altSeq is AGTG, returns 1 
    '''
    if len(refSeq) == len(altSeq):
        for index in xrange(len(refSeq)):
            if refSeq.lower()[index] != altSeq.lower()[index]:
                return index
            else:
                pass
    else:
        return "N/A"
            
def getZScore(maxEntScanScore, donor=False):
    '''
    Given a MaxEntScanScore, returns the zscore
    If donor is True, uses splice donor mean and std
    If donor is False, uses splice acceptor mean and std
    '''
    stdMeanData = json.load(open(os.path.join(os.path.dirname(__file__), 'brca.zscore.json')))
    if donor == False:
        std = stdMeanData["acceptors"]["std"]
        mean = stdMeanData["acceptors"]["mean"]
    else:
        std = stdMeanData["donors"]["std"]
        mean = stdMeanData["donors"]["mean"]
        
    zscore = (maxEntScanScore-mean)/std
    return zscore

def getRefAltScores(refSeq, altSeq, donor=False):
    '''
    Given ref and alt sequences and if sequence is in a splice donor region or not (True/False)
    Returns a dictionary containing raw MaxEntScan scores and zscores for ref and alt sequences
    '''
    if donor == False:
        refMaxEntScanScore = runMaxEntScan(refSeq, donor=False)
        refZScore = getZScore(refMaxEntScanScore, donor=False)
        altMaxEntScanScore = runMaxEntScan(altSeq, donor=False)
        altZScore = getZScore(altMaxEntScanScore, donor=False)    
    else:
        refMaxEntScanScore = runMaxEntScan(refSeq, donor=True)
        refZScore = getZScore(refMaxEntScanScore, donor=True)
        altMaxEntScanScore = runMaxEntScan(altSeq, donor=True)
        altZScore = getZScore(altMaxEntScanScore, donor=True)

    scoreDict = {"refScores": {"maxEntScanScore": refMaxEntScanScore,
                               "zScore": refZScore},
                 "altScores": {"maxEntScanScore": altMaxEntScanScore,
                               "zScore": altZScore}}
    return scoreDict

def getMaxEntScanScoresSlidingWindowSNS(variant, windowSize, donor=False):
    '''
    Given a variant and window size determines window sequences and scores for a sliding window
      that is the size of windowSize
    If donor=True, calculates MaxEntScan scores for splice donors
    If donor=False, calculates MaxEntScan scores for splice acceptors
    Returns a dictionary containing:
        1. window sequences - ref and alt seq for each window (variant in positions 1-windowSize)
        2. window scores - ref and alt MaxEntScan scores and zscores for each window
        3. window alt MaxEntScan scores - only contains alt MaxEntScan scores for each window
    '''
    varGenPos = int(variant["Pos"])
    varStrand = getVarStrand(variant)
    # use +- (windowSize - 1) to get (windowSize*2 - 1) bp region so that have sequence for:
    # each window of size windowSize bp with variant in each position (1-windowSize)
    # minus strand and plus strand are opposite for +- (windowSize - 1) to preserve sequence returned by getRefAltSeqs
    offset = windowSize - 1
    varPos = windowSize
    windowEnd = windowSize
    totalPositions = windowSize
    if varStrand == "-":
        regionStart = varGenPos + offset
        regionEnd = varGenPos - offset
    else:
        regionStart = varGenPos - offset
        regionEnd = varGenPos + offset
    refAltSeqs = getRefAltSeqs(variant, regionStart, regionEnd)
    refSeq = refAltSeqs["refSeq"]
    altSeq = refAltSeqs["altSeq"]
    windowStart = 0
    windowSeqs = {}
    windowScores = {}
    windowAltMaxEntScanScores = {}
    while windowStart < totalPositions:
        refWindowSeq = refSeq[windowStart:windowEnd]
        altWindowSeq = altSeq[windowStart:windowEnd]
        windowSeqs[varPos] = {"refSeq": refWindowSeq,
                              "altSeq": altWindowSeq}
        refAltWindowScores = getRefAltScores(refWindowSeq, altWindowSeq, donor=donor)
        windowScores[varPos] = {"refMaxEntScanScore": refAltWindowScores["refScores"]["maxEntScanScore"],
                                "refZScore": refAltWindowScores["refScores"]["zScore"],
                                "altMaxEntScanScore": refAltWindowScores["altScores"]["maxEntScanScore"],
                                "altZScore": refAltWindowScores["altScores"]["zScore"]}
        windowAltMaxEntScanScores[varPos] = refAltWindowScores["altScores"]["maxEntScanScore"]
        varPos -= 1
        windowStart += 1
        windowEnd += 1

    return {"windowSeqs": windowSeqs,
            "windowScores": windowScores,
            "windowAltMaxEntScanScores": windowAltMaxEntScanScores}

def getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, deNovoLength, donor=True, deNovo=False, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determines the maximum alt MaxEntScan score in 
       a sliding window of size STD_DONOR_SIZE with the variant in each position (1-STD_DONOR_SIZE) if donor = True
       a sliding window of size STD_ACC_SIZE with the variant in each position (1-STD_ACC_SIZE) if donor = False
    This function should be used to determine window in which de novo splicing is most likely to occur
    Function can only return highest scoring window details for either de novo donor OR de novo acceptor, not both
    If donor=True, function determines highest scoring window for potential de novo donor
    If donor=False, function determines highest scoring window for potential de novo acceptor
    Returns a dictionary containing the ref and alt MaxEntScan score and z-score and position of variant for the highest scoring window
    Ref and alt seqs for the highest scoring window are included in dictionary along with varStart (0-based index of variant for formatting)
       and varLength (equal to 1 for this function, becuase this function only works for single nucleotide substitution variants
    Dictionary also containing value "inExonicPortion" that has value either True or False
       If inExonicPortion = True, then variant is in length of bp specified by exonicPortionSize of highest scoring sliding window
       If inExonicPortion = False, then variant is NOT in length of bp specified by exonicPortionSize highest scoring sliding window 
    deNovoLength refers to the length of the exonic portion of a de novo splice acceptor
    deNovoDonorInRefAcc = False if NOT checking for de novo splice donors in reference splice acceptor sites
    deNovoDonorInRefAcc = True if checking for de novo splice donors in reference splice acceptor sites        
    '''
    if donor == True:
        # uses default window size for a splice donor region
        slidingWindowInfo = getMaxEntScanScoresSlidingWindowSNS(variant, STD_DONOR_SIZE, donor=donor)
    else:
        # uses default window size for a splice acceptor region
        slidingWindowInfo = getMaxEntScanScoresSlidingWindowSNS(variant, STD_ACC_SIZE, donor=donor)
    windowAltMaxEntScanScores = slidingWindowInfo["windowAltMaxEntScanScores"]
    # checks to see if variatn is within reference splice donor region
    inRefSpliceDonorRegion = varInSpliceRegion(variant, donor=True, deNovo=False)
    # checks to see if variant is within reference splice acceptor region
    inRefSpliceAccRegion = varInSpliceRegion(variant, donor=False, deNovo=False)
    # if variant in ref splice donor region (for de novo donor) or in ref splice acceptor region (for de novo acceptor),
    # then need to remove native splicing window from consideration for highest scoring window
    if (inRefSpliceDonorRegion == True or inRefSpliceAccRegion == True) and deNovoDonorInRefAcc == False:
        if donor == True:
            refSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            if getVarStrand(variant) == "+":
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["donorStart"], refSpliceBounds["donorEnd"], plusSeq=True)
            else:
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["donorStart"], refSpliceBounds["donorEnd"], plusSeq=False)
        else:
            refSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=True)
            deNovoOffset = deNovoLength - exonicPortionSize
            # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            if getVarStrand(variant) == "+":
                # acceptorEnd - deNovoOffset because genomic position increases from left to right on plus strand, refSeq reduced to correct length
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] - deNovoOffset), plusSeq=True)
            else:
                # acceptorEnd + deNovoOffset because genomic position decreases from left to right on minus strand, refSeq reduced to correct length
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] + deNovoOffset), plusSeq=False)
        for position, seqs in slidingWindowInfo["windowSeqs"].iteritems():
            if seqs["refSeq"] == refSpliceSeq:
                refSpliceWindow = position
                # removes reference splice window so it is not considered for de novo splicing
                del windowAltMaxEntScanScores[refSpliceWindow]
    # to get tuple containing sequence with variant position with maximum alt MaxEntScan score
    maxAltWindowScore = max(windowAltMaxEntScanScores.items(), key=lambda k: k[1])
    maxVarPosition = maxAltWindowScore[0]
    maxScores = slidingWindowInfo["windowScores"][maxVarPosition]
    maxSeqs = slidingWindowInfo["windowSeqs"][maxVarPosition]
    
    # determines if variant is in the exonic portion specified by exonicPortionLength
    inExonicPortion = False
    if donor == True:
    # determines if variant is in first exonicPortionSize bp of the donor region
        if maxVarPosition <= exonicPortionSize:
            inExonicPortion = True
    else:
    # determines if variant is in the last exonicPortionSize bp of the acceptor region
        if (STD_ACC_SIZE - maxVarPosition) < exonicPortionSize:
            inExonicPortion = True

    return {"refMaxEntScanScore": maxScores["refMaxEntScanScore"],
            "refZScore": maxScores["refZScore"],
            "altMaxEntScanScore": maxScores["altMaxEntScanScore"],
            "altZScore": maxScores["altZScore"],
            "refSeq": maxSeqs["refSeq"],
            "altSeq": maxSeqs["altSeq"],
            "varStart": maxVarPosition - 1,
            "varLength": 1,
            "varWindowPosition": maxVarPosition,
            "inExonicPortion": inExonicPortion}

def varInExonicPortion(variant, exonicPortionSize, deNovoLength, donor=True, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determines if variant in in the exonic portion as specified
    exonicPortionLength refers to the number of bases that are considered to be in the exon
    deNovoLength refers to the number of bases in the exon that are considered part of deNovo acceptor region
    if donor=True and exonicPortionSize=3, determines if variant is in first 3 bp of highest scoring window
    if donor=False and exonicPortionSize=3, determines if variant is in last 3 bp of highest scoring window
    If deNovoDonorInRefAcc=True, function is used in context of looking for de novo donor scores in ref splice acceptor sites
    if deNovoDonorInRefAcc=False, function is not used for de novo donors in ref splice acceptor sites
    Returns true if variant is in exonic portion, False otherwise
    '''
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize=exonicPortionSize,
                                                              deNovoLength=deNovoLength, donor=donor,
                                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    if slidingWindowInfo["inExonicPortion"] == True:
        return True
    return False

def getVarWindowPosition(variant, donor=True, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determines window position for highest scoring sliding window
    donor=True if function being used for splice donor, donor=False if function being used for splice acceptor
    Returns integer 1-STD_DONOR_SIZE based on variant position in highest scoring window if donor=True
    Returns integer 1-STD_ACC_SIZE based on variant position in highest scoring window if donor=False
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    '''
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                              donor=donor, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    varWindowPos = slidingWindowInfo["varWindowPosition"]
    return varWindowPos

def getClosestExonNumberIntronicSNS(variant, boundaries, donor=True):
    '''
    Given a variant and boundaries (either priors or enigma),
    1. Checks that variant is in an intron or UTR and is a SNS variant
    2. Determines the exon end that is closest to that variant
    Returns the closest exon end in the format "exonN"
    If variant is not in an intron or UTR, returns "exon0"
    '''
    varLoc = getVarLocation(variant, boundaries)
    if (varLoc == "intron_variant" or varLoc == "UTR_variant")  and getVarType(variant) == "substitution" and varInExon(variant) == False:
        exonBounds = getExonBoundaries(variant)
        varGenPos = variant["Pos"]
        exonIntronDiffs = {}
        for exon in exonBounds.keys():
            if getVarStrand(variant) == "+":
                if donor == True:
                    exonIntronDiff = int(varGenPos) - int(exonBounds[exon]["exonEnd"])
                else:
                    exonIntronDiff = int(exonBounds[exon]["exonStart"]) - int(varGenPos)
                if exonIntronDiff > 0:
                    exonIntronDiffs[exon] = exonIntronDiff
            else:
                if donor == True:
                    exonIntronDiff = int(exonBounds[exon]["exonEnd"]) - int(varGenPos)
                else:
                    exonIntronDiff = int(varGenPos) - int(exonBounds[exon]["exonStart"])
                if exonIntronDiff > 0:
                    exonIntronDiffs[exon] = exonIntronDiff
        closestExonInfo = min(exonIntronDiffs.items(), key=lambda k: k[1])
        return closestExonInfo[0]
    return "exon0"
                
def getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False, deNovoDonorInRefAcc=False, testMode=False):
    '''
    Given a variant, determines scores for closest reference splice sequence
    Also returns sequence of closest reference splice site
    deNovoOffset refers to difference between de novo acceptor length and exonic portion size
       If donor = True, looks for closest splice donor sequence
       If donor = False, looks for closest splice acceptor sequence
       If deNovo = True, accomodates for de novo splicing
         *Note only use argument deNovo=True in this function if donor=False
         *Function will not return correct sequence if donor=True and deNovo=True
    If exonic variant, returns a dictionary containing:
       MaxEntScan score, z-score, and splice site sequence for reference closest splice sequence
    If variant located in referene splice site, returns a dictionary containing:
       MaxEntScan score, z-score, and splice site sequence for that reference splice site sequence
    If intronic variant or variant in UTR, returns a dictionary containg:
       MaxEntScan score, z-score, and splice site sequence for reference closest splice site
       *Note if looking for closest ref acceptor for a variant in an intron, use deNovoOffset=0
    Return dictionary also contains necessary formatting variables for splice site sequence (exonStart, intronStart)
    deNovoDonorInRefAcc = False if NOT checking for de novo splice donor sites in reference splice acceptor sites
    deNovoDonorInRefAcc = True if checking for de novo splice donor sites in reference splice acceptor sites 
    '''
    varGenPos = int(variant["Pos"])
    varChrom = getVarChrom(variant)
    varLoc = getVarLocation(variant, "enigma")
    if (varInExon(variant) == True and deNovo == False) or (varLoc == "intron_variant" or varLoc == "UTR_variant"):
        if varInExon(variant) == True:
            exonNumber = getVarExonNumberSNS(variant)
            exonName = exonNumber
        if (varLoc == "intron_variant" or varLoc == "UTR_variant") and varInExon(variant) == False:
            exonNumber = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
            exonName = exonNumber
        if donor == True:
            refSpliceDonorBounds = getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
            closestSpliceBounds = refSpliceDonorBounds[exonNumber]
        else:
            refSpliceAccBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
            closestSpliceBounds = refSpliceAccBounds[exonNumber]
    if varInSpliceRegion(variant, donor=donor, deNovo=deNovo) == True and deNovoDonorInRefAcc == False:
        closestSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=deNovo)
        exonName = closestSpliceBounds["exonName"]
    if donor == True:
        if getVarStrand(variant) == "+":
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"], plusSeq=True)
        else:
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"], plusSeq=False)
        exonStart = 0
        intronStart = STD_EXONIC_PORTION
    else:
        # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
        if getVarStrand(variant) == "+":
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"], (closestSpliceBounds["acceptorEnd"] - deNovoOffset), plusSeq=True)
        else:
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"], (closestSpliceBounds["acceptorEnd"] + deNovoOffset), plusSeq=False)
        exonStart = len(refSeq) - STD_EXONIC_PORTION
        intronStart = 0
    if testMode == False:
        # to prevent issue with running max ent scan score on unittests
        if exonName == "exon0":
            return {"exonName": "N/A",
                    "sequence": "N/A",
                    "exonStart": "N/A",
                    "intronStart": "N/A",
                    "maxEntScanScore": "N/A",
                    "zScore": "N/A"}
        closestMaxEntScanScore = runMaxEntScan(refSeq, donor=donor)
        closestZScore = getZScore(closestMaxEntScanScore, donor=donor)
        return {"exonName": exonName,
                "sequence": refSeq.upper(),
                "exonStart": exonStart,
                "intronStart": intronStart,
                "maxEntScanScore": closestMaxEntScanScore,
                "zScore": closestZScore}
    else:
        return {"exonName": exonName,
                "sequence": refSeq.upper()}

def isCIDomainInRegion(regionStart, regionEnd, boundaries, gene):
    '''
    Given a region of interest, boundaries (either enigma or priors) and gene of interest
    Determines if there is an overlap between the region of interest and a CI domain
    For minus strand gene (BRCA1) regionStart > regionEnd
    For plus strand gene (BRCA2) regionStart < regionEnd
    Returns True if there is an overlap, False otherwise
    '''
    if gene == "BRCA1":
        for domain in brca1CIDomains[boundaries].keys():
            domainStart = brca1CIDomains[boundaries][domain]["domStart"]
            domainEnd = brca1CIDomains[boundaries][domain]["domEnd"]
            overlap = range(max(regionEnd, domainEnd), min(regionStart, domainStart) + 1)
            if len(overlap) > 0:
                return True
    elif gene == "BRCA2":
        for domain in brca2CIDomains[boundaries].keys():
            domainStart = brca2CIDomains[boundaries][domain]["domStart"]
            domainEnd = brca2CIDomains[boundaries][domain]["domEnd"]
            overlap = range(max(regionStart, domainStart), min(regionEnd, domainEnd) + 1)
            if len(overlap) > 0:
                return True
    return False

def getRefExonLength(variant, donor=True):
    '''
    Given a variant, returns the length of the reference exon
    If variant is in an exon, returns length of that exon
    If variant is in a reference splice region, returns length of exon in which exonic portion is included
    If variant is in intron, returns exon in which either closest splice donor or acceptor is included depending on donor argument
      If donor=True, returns exon length for previous exon
      If donor=False, returns exon length for subsequent exon
    '''
    if varInExon(variant) == True:
        varExonNum = getVarExonNumberSNS(variant)
    else:
        if varInSpliceRegion(variant, donor=donor, deNovo=False) == True:
            spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            varExonNum = spliceBounds["exonName"]
        else:
            varExonNum = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
    exonBounds = getExonBoundaries(variant)
    if getVarStrand(variant) == "-":
        varExonStart = int(exonBounds[varExonNum]["exonStart"])
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        exonLength = varExonStart - varExonEnd
    else:
        varExonStart = int(exonBounds[varExonNum]["exonStart"])
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        exonLength = varExonEnd - varExonStart
    return exonLength

def getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, exonicPortionSize, intronicPortionSize, donor=True):
    '''
    Given a variant's:
    genetic postion, strand, sliding window position with max MES score AND
      whether that position is within exonic portion of highest scoring window, exonic portion size, and intronic portion size
    Returns the position where splicing occurs for a de novo splice donor or acceptor (depending on donor=True argument)
    '''
    if varStrand == "+":
        if inExonicPortion == False:
            if donor == True:
                newSplicePos = int(varGenPos) - (varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) + abs(varWindowPos - intronicPortionSize)
        else:
            if donor == True:
                newSplicePos = int(varGenPos) + abs(varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) - (varWindowPos - intronicPortionSize)
    else:
        if inExonicPortion == False:
            if donor == True:
                newSplicePos = int(varGenPos) + (varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) - abs(varWindowPos - intronicPortionSize)
        else:
            if donor == True:
                newSplicePos = int(varGenPos) - abs(varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) + (varWindowPos - intronicPortionSize)
    return newSplicePos
    
def getAltExonLength(variant, exonicPortionSize, intronicPortionSize, deNovoDonorInRefAcc=False, donor=True):
    '''
    Given a variant and the exonic portion size,
    returns the length of the alternate exon after splicing occurs in max MES window
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    '''
    if varInExon(variant) == True:
        varExonNum = getVarExonNumberSNS(variant)
    else:
        if varInSpliceRegion(variant, donor=donor, deNovo=False) == True:
            spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            varExonNum = spliceBounds["exonName"]
        else:
            varExonNum = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
    exonBounds = getExonBoundaries(variant)
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH, donor=donor,
                                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    newSplicePos = getNewSplicePosition(variant["Pos"], getVarStrand(variant), slidingWindowInfo["varWindowPosition"],
                                        slidingWindowInfo["inExonicPortion"], exonicPortionSize, intronicPortionSize, donor=donor)
    if getVarStrand(variant) == "-":
        if donor == True:
            varExonStart = int(exonBounds[varExonNum]["exonStart"])
            # newSplicePos -1 to account for RefSeq numbering which starts to the right of the first base
            exonLength = varExonStart - (newSplicePos - 1)
        else:
            varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
            # varExonEnd +1 to account for RefSeq numbering which starts to the right of the first base
            exonLength = newSplicePos - (varExonEnd + 1)
    else:
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        if donor == True:
            refExonLength = getRefExonLength(variant, donor=donor)
            # need to compare to refExonLength because of + strand gene
            exonLength = refExonLength - (varExonEnd - newSplicePos)
        else:
            exonLength = varExonEnd - newSplicePos
    return exonLength

def compareRefAltExonLengths(refLength, altLength):
    '''
    Compares ref and alt exon lengths
    If both exon lengths % 3 are equal, then both ref and alt have the same reading frame
    Returns true if both ref and alt exon have the same reading frame, false otherwise
    '''
    if refLength % 3 == altLength % 3:
        return True
    else:
        return False

def isSplicingWindowInFrame(variant, exonicPortionSize, intronicPortionSize, deNovoDonorInRefAcc=False, donor=True):
    '''
    Given a variant, determines ref and alt exon length and compares them
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    intronicPortionSize referes to length in bp that is considered to be in intronic portion of splice site
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    If ref and alt exon are in the same reading frame, returns True
    '''
    refLength = getRefExonLength(variant, donor=donor)
    altLength = getAltExonLength(variant, exonicPortionSize, intronicPortionSize, deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=donor)
    return compareRefAltExonLengths(refLength, altLength)

def compareDeNovoWildTypeSplicePos(variant, exonicPortionSize, intronicPortionSize, deNovoDonorInRefAcc=False, donor=True):
    '''
    Given a variant, compares de novo splicing position with wild-type splicing position
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    intronicPortionSize referes to length in bp that is considered to be in intronic portion of splice site
    deNovoDonorInRefAcc argument=True if looking for de novo donor in reference splice acceptor region, False otherwise
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    If distance between de novo and wild-type donors is divisible by 3, returns True
    returns False otherwise
    '''
    if varInExon(variant) == True:
        varExonNum = getVarExonNumberSNS(variant)
    else:
        if varInSpliceRegion(variant, donor=donor, deNovo=False) == True:
            spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
            varExonNum = spliceBounds["exonName"]
        else:
            varExonNum = getClosestExonNumberIntronicSNS(variant, "enigma", donor=donor)
    varStrand = getVarStrand(variant)
    refExonBounds = getExonBoundaries(variant)
    slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH, donor=donor,
                                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    deNovoSplicePos = getNewSplicePosition(variant["Pos"], varStrand, slidingWindowInfo["varWindowPosition"],
                                           slidingWindowInfo["inExonicPortion"], exonicPortionSize, intronicPortionSize, donor=donor)
    if donor == True:
        wildTypeSplicePos = refExonBounds[varExonNum]["exonEnd"]
        if varStrand == "+":
            distanceBetween = wildTypeSplicePos - deNovoSplicePos
        else:
            wildTypeSplicePos = refExonBounds[varExonNum]["exonEnd"]
            # +1 because of RefSeq numbering
            distanceBetween = deNovoSplicePos - (wildTypeSplicePos + 1)
    else:
        wildTypeSplicePos = refExonBounds[varExonNum]["exonStart"]
        if varStrand == "+":
            distanceBetween = abs(deNovoSplicePos - wildTypeSplicePos)
        else:
            # +1 because of RefSeq numbering
            distanceBetween = abs((wildTypeSplicePos + 1) - deNovoSplicePos)
            
    if distanceBetween % 3 == 0:
        return True
    return False
        
def getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determines if there is a possibility of splice rescue
    deNovoDonorInRefAcc argument = True  if looking for deNovoDonor in ref acceptor site, False otherwise
    If there is a possibility of splice rescue, flags variant for further analysis
    Else assigns prior probability of pathogenecity and predicted qualitative ENIGMA class
    Extra flags are included to provide more information about why splicing rescue does or does not occur
      frameshiftFlag: equals 1 if variant causes a frameshift, 0 if not
      inExonicPortionFlag: equals 1 if variant IS in exonic portion of highest scoring window, 0 if not
      CIDomainRegionFlag: equals 1 if truncating region (between new donor and next reference splice acceptor) 
         includes a clinicall import domain, 0 if not
      isDivisibleFlag: equals 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3, 
         0 if not
    '''
    if varInIneligibleDeNovoExon(variant, donor=True) == True:
        return {"priorProb": 0.99,
                "enigmaClass": "class_5",
                "spliceRescue": 0,
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}
    # checks that variant causes a premature stop codon in an exon
    if getVarConsequences(variant) == "stop_gained" and varInExon(variant):
        spliceFlag = 0
        spliceRescue = 0
        frameshiftFlag = "-"
        inExonicPortionFlag = "-"
        CIDomainInRegionFlag = "-"
        isDivisibleFlag = "-"
        lowMESFlag = "-"
        # if variant is in specified exonic portion of highest scoring sliding window, no splice rescue
        if varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=True,
                              deNovoDonorInRefAcc=deNovoDonorInRefAcc) == True:
            priorProb = 0.97
            inExonicPortionFlag = 1
            spliceRescue = 0
        else:
            inFrame = isSplicingWindowInFrame(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=True)
            # if variant causes a frameshift, no splice rescue
            if inFrame == False:
                priorProb = 0.99
                frameshiftFlag = 1
                inExonicPortionFlag = 0
                spliceRescue = 0
            else:
                varGenPos = int(variant["Pos"])
                varStrand = getVarStrand(variant)
                varExonNum = getVarExonNumberSNS(variant)
                # varExonNum returns a string in the format "exonN"
                # nextExonNum parses out N from varExonNum and adds 1 to get next exon number key "exonN+1"
                # use [4:] to remove "exon" from "exonN" so can add 1 to N to get N+1
                nextExonNum = "exon" + str(int(varExonNum[4:]) + 1)
                # skips to exon 5 because exon 4 does not exist in BRCA1 refseq transcript
                if variant["Gene_Symbol"] == "BRCA1" and nextExonNum == "exon4":
                    nextExonNum = "exon5"
                refSpliceAccBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
                varWindowPos = getVarWindowPosition(variant, donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                inExonicPortion = varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=True,
                                                     deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                # gets region from new splice position to next splice acceptor
                regionStart = getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, STD_EXONIC_PORTION,
                                                   STD_ACC_INTRONIC_LENGTH, donor=True)
                regionEnd = refSpliceAccBounds[nextExonNum]["acceptorStart"]
                CIDomainInRegion = isCIDomainInRegion(regionStart, regionEnd, boundaries, variant["Gene_Symbol"])
                isDivisible = compareDeNovoWildTypeSplicePos(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                             deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=True)
                # if truncated region includes a clinically important domain or causes a frameshift
                if CIDomainInRegion == True:
                    priorProb = 0.99
                    spliceRescue = 0
                    CIDomainInRegionFlag = 1
                    inExonicPortionFlag = 0
                    frameshiftFlag = 0
                elif isDivisible == False:
                    priorProb = 0.99
                    spliceRescue = 0
                    isDivisibleFlag = 1
                    frameshiftFlag = 1
                    inExonicPortionFlag = 0
                    CIDomainInRegionFlag = 0
                else:
                    # possibility of splice rescue, check MES scores
                    frameshiftFlag = 0
                    inExonicPortionFlag = 0
                    CIDomainInRegionFlag = 0
                    isDivisibleFlag = 0
                    deNovoSpliceData = getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                                             donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                    altMES = deNovoSpliceData["altMaxEntScanScore"]
                    refZScore = deNovoSpliceData["refZScore"]
                    altZScore = deNovoSpliceData["altZScore"]
                    if altZScore <= refZScore:
                        # variant creates a weaker splice site than was there previously, no change in splicing
                        priorProb = 0.99
                        spliceRescue = 0
                        lowMESFlag = 1
                    else:
                        # variant creates a stronger splice site than reference sequence
                        if altMES < 6.2:
                            # still a weak splice site, so no change in splicing
                            priorProb = 0.99
                            spliceRescue = 0
                            lowMESFlag = 1
                        elif (altMES >= 6.2) and (altMES <= 8.5):
                            # still a weak splice site, but possibility of splice rescue
                            closestRefData = getClosestSpliceSiteScores(variant, boundaries, donor=True,
                                                                        deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                            closestZScore = closestRefData["zScore"]
                            if varInSpliceRegion(variant, donor=True, deNovo=False) == True:
                                closestAltData = getPriorProbRefSpliceDonorSNS(variant, boundaries)
                                closestZScore = closestAltData["altZScore"]
                            if altZScore > closestZScore:
                                # splice site created by variant is stronger than subsequent (closest) wild-type donor
                                priorProb = "N/A"
                                enigmaClass = "N/A"
                                spliceRescue = 1
                                spliceFlag = 1
                                lowMESFlag = 0
                            else:
                                # splice site created by variant is weaker than subsequent (closest) wild-type donor
                                priorProb = 0.99
                                spliceRescue = 0
                                lowMESFlag = 1
                        else:
                            # altMES > 8.5, strong splice site with higher possibility of splicing rescue
                            priorProb = "N/A"
                            enigmaClass = "N/A"
                            spliceRescue = 1
                            spliceFlag = 1
                            lowMESFlag = 0

        if priorProb == "N/A":
            pass
        else:
            enigmaClass = getEnigmaClass(priorProb)

        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag}

def getDeNovoSpliceFrameshiftStatus(variant, donor=True, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determiens if de novo splice site (either donor or acceptor based on donor argument)
      causes a frameshift
    Returns true if variant's de novo splice site causes a frameshift, false otherwise

    deNovoDonorInRefAcc is True if looking for de novo donor in reference splice acceptor site, False otherwise
    Frameshift is determined in 2 ways: inFrame and isDivisible
    '''
    # if inFrame == False then alt and ref exons are not in the same reading frame
    inFrame = isSplicingWindowInFrame(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                      deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=donor)
    # if isDivisble == Flase then distance between old and new splice position is not divislbe by 3
    isDivisible = compareDeNovoWildTypeSplicePos(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                 deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=donor)
    if inFrame == False or isDivisible == False:
        return True
    return False

def getEnigmaClass(priorProb):
    '''
    Given a prior probability of pathogenecity, returns a predicted qualitative ENIGMA class
    '''
    # if variant has prior prob = N/A then a predicted qualitative ENIGMA class will have already been determined
    if priorProb == "N/A":
        pass 
    else:
        if priorProb >= 0.99:
            return "class_5"
        elif priorProb <= 0.99 and priorProb >= 0.95:
            return "class_4"
        elif priorProb < 0.05 and priorProb >= 0.001:
            return "class_2"
        elif priorProb < 0.001:
            return "class_1"
        else:
            return "class_3"

def getPriorProbRefSpliceDonorSNS(variant, boundaries):
    '''
    Given a variant and location boundaries (either PRIORS or enigma)
    Checks that variant is in a splice donor site and is a single nucleotide substitution
    Returns a dictionary containing:
     prior probability of pathogenecity, predicted qualitative enigma class, ref and alt sequences, and ref and alt MES and zscores
     also contains variables needed to format ref and alt seq:
       - varStart: index (based on 0 indexing) where variant is located in altSeq
       - varLength: length of variant in altSeq (this function is only for SNS variants so varLength=1)
       - exonStart: denotes index (0-indexed) where exonic portion of sequence starts, always equal to 0 because splice donor is at end of exon
       - intronStart: denotes index (0-indexed) where intronic portion of sequence starts, equal to exonic portion size
     also has spliceSite variable = 1, so variant is marked as in a reference splice site
    '''
    varType = getVarType(variant)
    varLoc = getVarLocation(variant, boundaries)
    if varType == "substitution" and (varLoc == "splice_donor_variant" or varLoc == "CI_splice_donor_variant"):
        # to get region boundaries to get ref and alt seq
        spliceDonorBounds = getVarSpliceRegionBounds(variant, donor=True, deNovo=False)
        refAltSeqs = getRefAltSeqs(variant, spliceDonorBounds["donorStart"], spliceDonorBounds["donorEnd"])
        scores = getRefAltScores(refAltSeqs["refSeq"], refAltSeqs["altSeq"], donor=True)
        refMaxEntScanScore = scores["refScores"]["maxEntScanScore"]
        refZScore = scores["refScores"]["zScore"]
        altMaxEntScanScore = scores["altScores"]["maxEntScanScore"]
        altZScore = scores["altScores"]["zScore"]
        if altMaxEntScanScore >= refMaxEntScanScore:
            priorProb = 0.04
        elif (refZScore < -1.5) and ((refZScore - altZScore) > 0.5):
            priorProb = 0.97
        elif (refZScore < -1.5) and ((refZScore - altZScore) < 0.5):
            priorProb = 0.34
        else:
            if altZScore > 0.0:
                priorProb = 0.04
            elif altZScore <= 0.0 and altZScore >= -2:
                priorProb = 0.34
            else:
                priorProb = 0.97    
        enigmaClass = getEnigmaClass(priorProb)
        
        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "refSeq": refAltSeqs["refSeq"].upper(),
                "altSeq": refAltSeqs["altSeq"].upper(),
                "varStart": getVarSeqIndexSNS(refAltSeqs["refSeq"], refAltSeqs["altSeq"]),
                "varLength": 1,
                "exonStart": 0,
                "intronStart": STD_EXONIC_PORTION,
                "refMaxEntScanScore": refMaxEntScanScore,
                "altMaxEntScanScore": altMaxEntScanScore,
                "refZScore": refZScore,
                "altZScore": altZScore,
                "spliceSite": 1}
    
def getPriorProbRefSpliceAcceptorSNS(variant, boundaries):
    '''
    Given a variant and location boundaries (either PRIORS or enigma)
    Checks that variant is in a splice acceptor site and is a single nucleotide substitution
    Returns a dictionary containing:
     prior probability of pathogenecity, predicted qualitative enigma class, ref and alt sequences, and ref and alt MES and zscores
     also contains variables needed to format ref and alt seq:
       - varStart: index (based on 0 indexing) where variant is located in altSeq
       - varLength: length of variant in altSeq (this function is only for SNS variants so varLength=1)
       - exonStart: denotes index (0-indexed) where exonic portion of sequence starts, always equal to 0 because splice donor is at end of exon
       - intronStart: denotes index (0-indexed) where intronic portion of sequence starts, equal to exonic portion size
     also has spliceSite variable = 1, so variant is marked as in a reference splice site
    '''
    varType = getVarType(variant)
    varLoc = getVarLocation(variant, boundaries)
    if varType == "substitution" and (varLoc == "splice_acceptor_variant" or varLoc == "CI_splice_acceptor_variant"):
        # to get region boundaires to get ref and alt seq
        spliceAcceptorBounds = getVarSpliceRegionBounds(variant, donor=False, deNovo=False)
        refAltSeqs = getRefAltSeqs(variant, spliceAcceptorBounds["acceptorStart"], spliceAcceptorBounds["acceptorEnd"])
        scores = getRefAltScores(refAltSeqs["refSeq"], refAltSeqs["altSeq"], donor=False)
        refMaxEntScanScore = scores["refScores"]["maxEntScanScore"]
        refZScore = scores["refScores"]["zScore"]
        altMaxEntScanScore = scores["altScores"]["maxEntScanScore"]
        altZScore = scores["altScores"]["zScore"]
        if altMaxEntScanScore >= refMaxEntScanScore:
            priorProb = 0.04
        elif (refZScore < -1.0) and ((refZScore - altZScore) > 0.5):
            priorProb = 0.97
        elif (refZScore < -1.0) and ((refZScore - altZScore) < 0.5):
            priorProb = 0.34
        else:
            if altZScore > 0.5:
                priorProb = 0.04
            elif altZScore <= 0.5 and altZScore >= -1.5:
                priorProb = 0.34
            else:
                priorProb = 0.97
        enigmaClass = getEnigmaClass(priorProb)

        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "refSeq": refAltSeqs["refSeq"].upper(),
                "altSeq": refAltSeqs["altSeq"].upper(),
                "varStart": getVarSeqIndexSNS(refAltSeqs["refSeq"], refAltSeqs["altSeq"]),
                "varLength": 1,
                "exonStart": len(refAltSeqs["refSeq"]) - STD_EXONIC_PORTION,
                "intronStart": 0,
                "refMaxEntScanScore": refMaxEntScanScore,
                "altMaxEntScanScore": altMaxEntScanScore,
                "refZScore": refZScore,
                "altZScore": altZScore,
                "spliceSite": 1}

def getPriorProbAfterGreyZoneSNS(variant, boundaries):
    '''
    Given a variant and location boundaries (either PRIORS or enigma)
    Checks that variant is after the grey zone and is a single nucleotide substitution
    Checks that variant is either a missense or nonsense mutation
    Returns a dictionary containing prior probability of pathogenecity and predicted qualitative enigma class
    Dictionary contains other values that are set to either N/A, "-", or 0 because they are not relevant
    '''
    varType = getVarType(variant)
    varLoc = getVarLocation(variant, boundaries)
    if varType == "substitution" and varLoc == "after_grey_zone_variant":
        varCons = getVarConsequences(variant)
        if varCons == "stop_gained" or varCons == "missense_variant" or varCons == "synonymous_variant":
            priorProb = "N/A"
            enigmaClass = "class_2"
        return {"applicablePrior": priorProb,
                "applicableEnigmaClass": enigmaClass,
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
                "refRefDonorMES": "-",
                "refRefDonorZ": "-",
                "altRefDonorMES": "-",
                "altRefDonorZ": "-",
                "refRefDonorSeq": "-",
                "altRefDonorSeq": "-",
                "refDonorVarStart": "-",
                "refDonorVarLength": "-",
                "refDonorExonStart": "-",
                "refDonorIntronStart": "-",
                "refDeNovoDonorMES": "-",
                "refDeNovoDonorZ": "-",
                "altDeNovoDonorMES": "-",
                "altDeNovoDonorZ": "-",
                "refDeNovoDonorSeq": "-",
                "altDeNovoDonorSeq": "-",
                "deNovoDonorVarStart": "-",
                "deNovoDonorVarLength": "-",
                "deNovoDonorExonStart": "-",
                "deNovoDonorIntronStart": "-",
                "closestDonorRefMES": "-",
                "closestDonorRefZ": "-",
                "closestDonorRefSeq": "-",
                "closestDonorAltMES": "-",
                "closestDonorAltZ": "-",
                "closestDonorAltSeq": "-",
                "closestDonorExonStart": "-",
                "closestDonorIntronStart": "-",
                "deNovoDonorAltGreaterRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestAltFlag": "N/A",
                "deNovoDonorFrameshiftFlag": "N/A",
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
                "refRefAccMES": "-",
                "refRefAccZ": "-",
                "altRefAccMES": "-",
                "altRefAccZ": "-",
                "refRefAccSeq": "-",
                "altRefAccSeq": "-",
                "refAccVarStart": "-",
                "refAccVarLength": "-",
                "refAccExonStart": "-",
                "refAccIntronStart": "-",                
                "refDeNovoAccMES": "-",
                "refDeNovoAccZ": "-",
                "altDeNovoAccMES": "-",
                "altDeNovoAccZ": "-",
                "refDeNovoAccSeq": "-",
                "altDeNovoAccSeq": "-",
                "deNovoAccVarStart": "-",
                "deNovoAccVarLength": "-",
                "deNovoAccExonStart": "-",
                "deNovoAccIntronStart": "-",
                "closestAccRefMES": "-",
                "closestAccRefZ": "-",
                "closestAccRefSeq": "-",
                "closestAccAltMES": "-",
                "closestAccAltZ": "-",
                "closestAccAltSeq": "-",
                "closestAccExonStart": "-",
                "closestAccIntronStart": "-",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}

def varInIneligibleDeNovoExon(variant, donor=True):
    '''
    Given a variant and donor argument:
        (donor=True, if for a de novo donor, donor=False, if for a de novo acceptor)
    Determines whether that variant is in an eligible exon to be evaluated for de novo splicing
    If in ineligible exon, returns True, returns False otherwise
    '''
    if varInExon(variant) == True:
        varGene = variant["Gene_Symbol"]
        varExon = getVarExonNumberSNS(variant)
        if donor == True:
            # last exon not eligible for de novo splice donor
            if varGene == "BRCA1" and varExon == "exon24":
                return True
            elif varGene == "BRCA2" and varExon == "exon27":
                return True
        else:
            # first exon not eligible for de novo splice acceptor
            if varExon == "exon1":
                return True
        return False

def getPriorProbDeNovoDonorSNS(variant, exonicPortionSize, deNovoDonorInRefAcc=False):
    '''
    Given a variant and exonicPortionSize
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in an exon or is in a reference splice donor region
    Returns a dictionary containing: 
      prior probability of pathogenecity and predicted qualitative engima class 
      deNovo donor MaxEntScan scores, zscores, and sequences for ref and alt
      closest donor MaxEntScan score, zscore, and sequence
      altGreaterRefFlag: equals 1 if altZScore > refZScore, 0 otherwise
      altGreaterClosestRefFlag: equals 1 if altZScore > closestRefZScore, 0 otherwise
      altGreaterClosestAltFlag: equals 1 if variant is in a native donor and altZScore > closestAltZScore, 0 if not,
                                N/A if not in a native donor
      frameshiftFlag: equals 1 if new splice site causes a frameshift
      variables necessary for formatting any sequences that are returned by this function
    deNovoDonorInRefAcc = False if NOT checking for de novo donor in reference splice acceptor site
    deNovoDonorInRefAcc = True if checking for de novo donors in reference splice acceptor site
    '''
    if getVarType(variant) == "substitution":
        if varInSpliceRegion(variant, donor=True, deNovo=True) == True:
            if varInExon(variant) == True and varInIneligibleDeNovoExon(variant, donor=True) == True:
                return {"priorProb": "N/A",
                        "enigmaClass": "N/A",
                        "refMaxEntScanScore": "N/A",
                        "altMaxEntScanScore": "N/A",
                        "refZScore": "N/A",
                        "altZScore": "N/A",
                        "refSeq": "N/A",
                        "altSeq": "N/A",
                        "varStart": "N/A",
                        "varLength": "N/A",
                        "exonStart": "N/A",
                        "intronStart": "N/A",
                        "closestRefMaxEntScanScore": "N/A",
                        "closestRefZScore": "N/A",
                        "closestRefSeq": "N/A",
                        "closestAltMaxEntScanScore": "N/A",
                        "closestAltZScore": "N/A",
                        "closestAltSeq": "N/A",
                        "closestExonStart": "N/A",
                        "closestIntronStart": "N/A",
                        "frameshiftFlag": "N/A",
                        "altGreaterRefFlag": "N/A",
                        "altGreaterClosestRefFlag": "N/A",
                        "altGreaterClosestAltFlag": "N/A"}
            slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH,
                                                                      donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            subDonorInfo = getClosestSpliceSiteScores(variant, STD_DE_NOVO_OFFSET, donor=True, deNovo=False,
                                                      deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            subZScore = subDonorInfo["zScore"]
            refAltZScore = "N/A"
            refAltMES = "N/A"
            refAltSeq = "N/A"
            altGreaterClosestRefFlag = 0
            altGreaterClosestAltFlag = "N/A"
            if varInSpliceRegion(variant, donor=True, deNovo=False) == True:
                refDonorInfo = getPriorProbRefSpliceDonorSNS(variant, "enigma")
                refAltZScore = refDonorInfo["altZScore"]
                refAltMES = refDonorInfo["altMaxEntScanScore"]
                refAltSeq = refDonorInfo["altSeq"]
                altGreaterClosestAltFlag = 0
            frameshiftFlag = 0
            frameshiftStatus = getDeNovoSpliceFrameshiftStatus(variant, donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            if frameshiftStatus == True:
                frameshiftFlag = 1
            altZScore = slidingWindowInfo["altZScore"]
            refZScore = slidingWindowInfo["refZScore"]
            altGreaterRefFlag = 0
            if altZScore <= refZScore:
                priorProb = 0
            else:
                altGreaterRefFlag = 1
                if altZScore < -2.0:
                    priorProb = 0.02
                elif altZScore >= -2.0 and altZScore < 0.0:
                    priorProb = 0.3
                else:
                    priorProb = 0.64
            if (altZScore > subZScore and refAltZScore == "N/A") or (altZScore > refAltZScore and refAltZScore != "N/A"):
                # promote prior prob by one step
                if priorProb == 0:
                    priorProb = 0.3
                elif priorProb == 0.02:
                    priorProb = 0.3
                elif priorProb == 0.3:
                    priorProb = 0.64
                else:
                    priorProb = priorProb

            if altZScore > subZScore:
                altGreaterClosestRefFlag = 1
            if altZScore > refAltZScore and refAltZScore != "N/A":
                altGreaterClosestAltFlag = 1
            
            if priorProb == 0: 
                priorProb = "N/A"
                enigmaClass = "N/A"
            else:
                priorProb = priorProb
                enigmaClass = getEnigmaClass(priorProb)
                
            return {"priorProb": priorProb,
                    "enigmaClass": enigmaClass,
                    "refMaxEntScanScore": slidingWindowInfo["refMaxEntScanScore"],
                    "altMaxEntScanScore": slidingWindowInfo["altMaxEntScanScore"],
                    "refZScore": refZScore,
                    "altZScore": altZScore,
                    "refSeq": slidingWindowInfo["refSeq"],
                    "altSeq": slidingWindowInfo["altSeq"],
                    "varStart": slidingWindowInfo["varStart"],
                    "varLength": slidingWindowInfo["varLength"],
                    "exonStart": 0,
                    "intronStart": STD_EXONIC_PORTION,
                    "closestRefMaxEntScanScore": subDonorInfo["maxEntScanScore"],
                    "closestRefZScore": subDonorInfo["zScore"],
                    "closestRefSeq": subDonorInfo["sequence"],
                    "closestAltMaxEntScanScore": refAltMES,
                    "closestAltZScore": refAltZScore,
                    "closestAltSeq": refAltSeq,
                    "closestExonStart": subDonorInfo["exonStart"],
                    "closestIntronStart": subDonorInfo["intronStart"],
                    "altGreaterRefFlag": altGreaterRefFlag,
                    "altGreaterClosestRefFlag": altGreaterClosestRefFlag,
                    "altGreaterClosestAltFlag": altGreaterClosestAltFlag,
                    "frameshiftFlag": frameshiftFlag}

def getPriorProbDeNovoAcceptorSNS(variant, exonicPortionSize, deNovoLength):
    '''
    Given a variant, exonic portion size, and de novo length:
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in de novo splice acceptor region 
         de novo splice acceptor region defined by deNovoLength
    Returns a dictionary containing: 
      prior probability of pathogenecity and predicted qualitative engima class (both N/A)
      deNovo acceptor MaxEntScan scores, zscores, and sequence for ref and alt
      MaxEntScan score, zscore, and sequence for closest ref splice acceptor
      frameshiftFlag: equals 1 if new splice site causes a frameshift
      altGreaterRefFlag: which equals 1 altZScore > refZScore, 0 otherwise
      altGreaterClosestRefFlag: which equals 1 if altZScore > closestRefZScore, 0 otherwise
      altGreaterClosestAltFlag: which equals 1 if variant in native acceptor and altZScore > closestAltZScore, if not,
                                N/A if not in a native acceptor
      variables necessary for formatting the returned sequences
    '''
    if getVarType(variant) == "substitution":
        if varInSpliceRegion(variant, donor=False, deNovo=True) == True:
            slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, deNovoLength,
                                                                        donor=False)
            deNovoOffset = deNovoLength - exonicPortionSize
            closestAccInfo = getClosestSpliceSiteScores(variant, deNovoOffset, donor=False, deNovo=True)
            refAltZScore = "N/A"
            refAltMES = "N/A"
            refAltSeq = "N/A"
            altGreaterClosestRefFlag = 0
            altGreaterClosestAltFlag = "N/A"
            if varInSpliceRegion(variant, donor=False, deNovo=False) == True:
                refAccInfo = getPriorProbRefSpliceAcceptorSNS(variant, "enigma")
                refAltZScore = refAccInfo["altZScore"]
                refAltMES = refAccInfo["altMaxEntScanScore"]
                refAltSeq = refAccInfo["altSeq"]
                altGreaterClosestAltFlag = 0
            frameshiftFlag = 0
            frameshiftStatus = getDeNovoSpliceFrameshiftStatus(variant, donor=False, deNovoDonorInRefAcc=False)
            if frameshiftStatus == True:
                frameshiftFlag = 1
            altZScore = slidingWindowInfo["altZScore"]
            refZScore = slidingWindowInfo["refZScore"]
            closestZScore = closestAccInfo["zScore"]
            altGreaterRefFlag = 0
            if altZScore > refZScore:
                altGreaterRefFlag = 1
            if altZScore > closestZScore:
                altGreaterClosestRefFlag = 1
            if altZScore > refAltZScore and refAltZScore != "N/A":
                altGreaterClosestAltFlag = 1
            return {"priorProb": "N/A",
                    "enigmaClass": "N/A",
                    "refMaxEntScanScore": slidingWindowInfo["refMaxEntScanScore"],
                    "altMaxEntScanScore": slidingWindowInfo["altMaxEntScanScore"],
                    "refZScore": refZScore,
                    "altZScore": altZScore,
                    "refSeq": slidingWindowInfo["refSeq"],
                    "altSeq": slidingWindowInfo["altSeq"],
                    "varStart": slidingWindowInfo["varStart"],
                    "varLength": slidingWindowInfo["varLength"],
                    "exonStart": STD_ACC_SIZE - STD_EXONIC_PORTION,
                    "intronStart": 0,
                    "closestRefMaxEntScanScore": closestAccInfo["maxEntScanScore"],
                    "closestRefZScore": closestAccInfo["zScore"],
                    "closestRefSeq": closestAccInfo["sequence"],
                    "closestAltMaxEntScanScore": refAltMES,
                    "closestAltZScore": refAltZScore,
                    "closestAltSeq": refAltSeq,
                    "closestExonStart": closestAccInfo["exonStart"],
                    "closestIntronStart": closestAccInfo["intronStart"],
                    "altGreaterRefFlag": altGreaterRefFlag,
                    "altGreaterClosestRefFlag": altGreaterClosestRefFlag,
                    "altGreaterClosestAltFlag": altGreaterClosestAltFlag,
                    "frameshiftFlag": frameshiftFlag}

def getPriorProbSpliceDonorSNS(variant, boundaries, variantData):
    '''
    Given a variant, boundaries (either PRIORS or ENIGMA), and a list of dictionaries with variant data
    Determines reference donor and de novo donor scores for variant
    If variant causes a nonsense mutation, determines if splice rescue occurs
    Returns dicitionary containing scores and sequences for ref and de novo splice donor/acceptor,
        closest splice site scores and sequences, and protein prior if variant in exon
        "N/A" if score, sequence or flag not applicable for variant
    Also contains other values:
        applicable prior, highest prior if variant has multiple priors
        ref prior, prior prob for reference splice donor
        de novo prior, prior prob for de novo donor sequence
        splice flag = 1, because variant in reference splice site
        deNovoDonorFrameshiftFlag = 1 if new splice site causes a frameshift
        deNovoDonorAltGreaterRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestAltFlag = 1 if variant is possible de novo donor variant
        deNovoAccFrameshiftFlag = N/A, because not applicable for variants in ref donor sites
        deNovoAccAltGreaterRefFlag = N/A, because not applicable for variants in ref donor sites
        deNovoAccAltGreaterClosestRefFlag = N/A, because not applicable for variants in ref donor sites
        deNovoAccAltGreaterClosestAltFlag = N/A, because not applicable for variants in ref donor sites
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise, N/A if not nonsense variant
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshiftFlag = 1 if nonsense variant causes a frameshift mutation also, 0 if not, N/A if not nonsense variant
        inExonicPortionFlag = 1 if variant IS in exonic portion of highest scoring window and variant is a nonsense variant, 
           0 if NOT in exonic portion, N/A if not nonsense variant
        CIDomainRegionFlag = 1 if truncating region (between new donor and next reference splice acceptor) 
           includes a clinically import domain for nonsense variant, 0 if not, N/A if not nonsense variant
        isDivisibleFlag = 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3 for nonsense variant, 
           0 if it IS divisible, N/A if not a nonsense variant
    Dictionary also contains formatting variables for each listed sequence
    ''' 
    if varInSpliceRegion(variant, donor=True, deNovo=False) and getVarType(variant) == "substitution":
        refSpliceInfo = getPriorProbRefSpliceDonorSNS(variant, boundaries)
        deNovoSpliceInfo = getPriorProbDeNovoDonorSNS(variant, STD_EXONIC_PORTION)
        deNovoPrior = deNovoSpliceInfo["priorProb"]
        refPrior = refSpliceInfo["priorProb"]
        proteinPrior = "N/A"
        if varInExon(variant) == True:
            proteinInfo = getPriorProbProteinSNS(variant, variantData)
            proteinPrior = proteinInfo["priorProb"]
        if deNovoPrior != "N/A" and proteinPrior != "N/A":
            applicablePrior = max(deNovoPrior, refPrior, proteinPrior)
        elif deNovoPrior == "N/A" and proteinPrior != "N/A":
            applicablePrior = max(refPrior, proteinPrior)
        elif deNovoPrior != "N/A" and proteinPrior == "N/A":
            applicablePrior = max(deNovoPrior, refPrior)
        elif deNovoPrior == "N/A" and proteinPrior == "N/A":
            applicablePrior = refPrior

        spliceRescue = "N/A"
        spliceFlag = 0
        frameshiftFlag = "N/A"
        inExonicPortionFlag = "N/A"
        CIDomainInRegionFlag = "N/A"
        isDivisibleFlag = "N/A"
        lowMESFlag = "N/A"
        # to check for nonsense variants in exonic portion of splice donor site
        if varInExon(variant) == True and getVarConsequences(variant) == "stop_gained":
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries)
            applicablePrior = nonsenseData["priorProb"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshiftFlag = nonsenseData["frameshiftFlag"]
            inExonicPortionFlag = nonsenseData["inExonicPortionFlag"]
            CIDomainInRegionFlag = nonsenseData["CIDomainInRegionFlag"]
            isDivisibleFlag = nonsenseData["isDivisibleFlag"]
            lowMESFlag = nonsenseData["lowMESFlag"]
            
        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": getEnigmaClass(applicablePrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": refPrior,
                "deNovoDonorPrior": deNovoPrior,
                "refRefDonorMES": refSpliceInfo["refMaxEntScanScore"],
                "refRefDonorZ": refSpliceInfo["refZScore"],
                "altRefDonorMES": refSpliceInfo["altMaxEntScanScore"],
                "altRefDonorZ": refSpliceInfo["altZScore"],
                "refRefDonorSeq": refSpliceInfo["refSeq"],
                "altRefDonorSeq": refSpliceInfo["altSeq"],
                "refDonorVarStart": refSpliceInfo["varStart"],
                "refDonorVarLength": refSpliceInfo["varLength"],
                "refDonorExonStart": refSpliceInfo["exonStart"],
                "refDonorIntronStart": refSpliceInfo["intronStart"],
                "refDeNovoDonorMES": deNovoSpliceInfo["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoSpliceInfo["refZScore"],
                "altDeNovoDonorMES": deNovoSpliceInfo["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoSpliceInfo["altZScore"],
                "refDeNovoDonorSeq": deNovoSpliceInfo["refSeq"],
                "altDeNovoDonorSeq": deNovoSpliceInfo["altSeq"],
                "deNovoDonorVarStart": deNovoSpliceInfo["varStart"],
                "deNovoDonorVarLength": deNovoSpliceInfo["varLength"],
                "deNovoDonorExonStart": deNovoSpliceInfo["exonStart"],
                "deNovoDonorIntronStart": deNovoSpliceInfo["intronStart"],
                "closestDonorRefMES": deNovoSpliceInfo["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoSpliceInfo["closestRefZScore"],
                "closestDonorRefSeq": deNovoSpliceInfo["closestRefSeq"],
                "closestDonorAltMES": deNovoSpliceInfo["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoSpliceInfo["closestAltZScore"],
                "closestDonorAltSeq": deNovoSpliceInfo["closestAltSeq"],
                "closestDonorExonStart": deNovoSpliceInfo["closestExonStart"],
                "closestDonorIntronStart": deNovoSpliceInfo["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoSpliceInfo["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoSpliceInfo["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoSpliceInfo["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoSpliceInfo["frameshiftFlag"],
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": "N/A",
                "refDeNovoAccZ": "N/A",
                "altDeNovoAccMES": "N/A",
                "altDeNovoAccZ": "N/A",
                "refDeNovoAccSeq": "N/A",
                "altDeNovoAccSeq": "N/A",
                "deNovoAccVarStart": "N/A",
                "deNovoAccVarLength": "N/A",
                "deNovoAccExonStart": "N/A",
                "deNovoAccIntronStart": "N/A",
                "closestAccRefMES": "N/A",
                "closestAccRefZ": "N/A",
                "closestAccRefSeq": "N/A",
                "closestAccAltMES": "N/A",
                "closestAccAltZ": "N/A",
                "closestAccAltSeq": "N/A",
                "closestAccExonStart": "N/A",
                "closestAccIntronStart": "N/A",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": refSpliceInfo["spliceSite"],
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag}

def getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData):
    '''
    Given a variant, boundaries (either PRIORS or ENIGMA), and list of dictionaries with variant data
    Determines reference and de novo acceptor scores for variant
      If variant in exon, also determines de novo donor scores and protein prior
    If variant causes a nonsense mutation, determines if splice rescue occurs
    Returns dicitionary containing scores for and sequences ref and de novo splice donor/acceptor,
        closest splice site scores and sequences and protein prior if variant in exon
        "N/A" if score, sequence, or flag not applicable for variant
    Also contains other values:
        applicable prior, highest prior if variant has multiple priors
        applicable classe, highest predicted qualitative enigma class based on applicable prior
        ref prior, prior prob for reference splice sequence
        de novo donor and acceptor priors, prior prob for de novo splice sequence
        splice flag = 1, because variant in reference splice site
        deNovoAccFrameshiftFlag = 1 if de novo acceptor causes a frameshift
        deNovoAccAltGreaterRefFlag = 1 if variant is possible de novo acceptor variant
        deNovoAccAltGreaterClosestRefFlag = 1 if variant is possible de novo acceptor variant
        deNovoAccAltGreaterClosestAltFlag = 1 if variant is possible de novo acceptor variant
        deNovoDonroFrameshiftFlag =1 if de novo donor causes a frameshift
        deNovoDonorAltGreaterRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestRefFlag = 1 if variant is possible de novo donor variant
        deNovoDonorAltGreaterClosestAltFlag = 1 if variant is possible de novo donor variant
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise, N/A if not nonsense variant
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshiftFlag = 1 if nonsense variant causes a frameshift mutation also, 0 if not, N/A if not nonsense variant
        inExonicPortionFlag = 1 if variant IS in exonic portion of highest scoring window and variant is a nonsense variant, 
           0 if NOT in exonic portion, N/A if not nonsense variant
        CIDomainRegionFlag = 1 if truncating region (between new donor and next reference splice acceptor) 
           includes a clinically import domain for nonsense variant, 0 if not, N/A if not nonsense variant
        isDivisibleFlag = 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3 for nonsense variant, 
           0 if it IS divisible, N/A if not a nonsense variant
    Dictionary also contains formatting variables for each listed sequence
    '''
    if varInSpliceRegion(variant, donor=False, deNovo=False) and getVarType(variant) == "substitution":
        refSpliceInfo = getPriorProbRefSpliceAcceptorSNS(variant, boundaries)
        deNovoAccInfo = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        refPrior = refSpliceInfo["priorProb"]
        proteinPrior = "N/A"
        applicablePrior = refSpliceInfo["priorProb"]
        if varInExon(variant) == True:
            deNovoDonorInfo = getPriorProbDeNovoDonorSNS(variant, STD_EXONIC_PORTION, deNovoDonorInRefAcc=True)
            deNovoDonorPrior = deNovoDonorInfo["priorProb"]
            proteinInfo = getPriorProbProteinSNS(variant, variantData)
            proteinPrior = proteinInfo["priorProb"]
            if deNovoDonorPrior != "N/A" and proteinPrior != "N/A":
                applicablePrior = max(deNovoDonorPrior, proteinPrior, refPrior)
            else:
                applicablePrior = max(proteinPrior, refPrior)
        else:
            deNovoDonorPrior = "N/A"
            deNovoDonorInfo = {"refMaxEntScanScore": "N/A",
                               "refZScore": "N/A",
                               "altMaxEntScanScore": "N/A",
                               "altZScore": "N/A",
                               "refSeq": "N/A",
                               "altSeq": "N/A",
                               "varStart": "N/A",
                               "varLength": "N/A",
                               "exonStart": "N/A",
                               "intronStart": "N/A",
                               "closestRefMaxEntScanScore": "N/A",
                               "closestRefZScore": "N/A",
                               "closestRefSeq": "N/A",
                               "closestAltMaxEntScanScore": "N/A",
                               "closestAltZScore": "N/A",
                               "closestAltSeq": "N/A",
                               "closestExonStart": "N/A",
                               "closestIntronStart": "N/A",
                               "altGreaterRefFlag": "N/A",
                               "altGreaterClosestRefFlag": "N/A",
                               "altGreaterClosestAltFlag": "N/A",
                               "frameshiftFlag": "N/A",
                               "priorProb": "N/A",
                               "enigmaClass": "N/A"}

        spliceRescue = "N/A"
        spliceFlag = 0
        frameshiftFlag = "N/A"
        inExonicPortionFlag = "N/A"
        CIDomainInRegionFlag = "N/A"
        isDivisibleFlag = "N/A"
        lowMESFlag = "N/A"
        # to check for nonsense variants in exonic portion of splice acceptor site
        if varInExon(variant) == True and getVarConsequences(variant) == "stop_gained":
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=True)
            applicablePrior = nonsenseData["priorProb"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshiftFlag = nonsenseData["frameshiftFlag"]
            inExonicPortionFlag = nonsenseData["inExonicPortionFlag"]
            CIDomainInRegionFlag = nonsenseData["CIDomainInRegionFlag"]
            isDivisibleFlag = nonsenseData["isDivisibleFlag"]
            lowMESFlag = nonsenseData["lowMESFlag"]

        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": getEnigmaClass(applicablePrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": deNovoDonorPrior,
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorInfo["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorInfo["refZScore"],
                "altDeNovoDonorMES": deNovoDonorInfo["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorInfo["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorInfo["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorInfo["altSeq"],
                "deNovoDonorVarStart": deNovoDonorInfo["varStart"],
                "deNovoDonorVarLength": deNovoDonorInfo["varLength"],
                "deNovoDonorExonStart": deNovoDonorInfo["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorInfo["intronStart"],
                "closestDonorRefMES": deNovoDonorInfo["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorInfo["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorInfo["closestRefSeq"],
                "closestDonorAltMES": deNovoDonorInfo["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoDonorInfo["closestAltZScore"],
                "closestDonorAltSeq": deNovoDonorInfo["closestAltSeq"],
                "closestDonorExonStart": deNovoDonorInfo["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorInfo["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorInfo["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorInfo["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorInfo["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorInfo["frameshiftFlag"],
                "deNovoAccPrior": deNovoAccInfo["priorProb"],
                "refAccPrior": refPrior,
                "refRefAccMES": refSpliceInfo["refMaxEntScanScore"],
                "refRefAccZ": refSpliceInfo["refZScore"],
                "altRefAccMES": refSpliceInfo["altMaxEntScanScore"],
                "altRefAccZ": refSpliceInfo["altZScore"],
                "refRefAccSeq": refSpliceInfo["refSeq"],
                "altRefAccSeq": refSpliceInfo["altSeq"],
                "refAccVarStart": refSpliceInfo["varStart"],
                "refAccVarLength": refSpliceInfo["varLength"],
                "refAccExonStart": refSpliceInfo["exonStart"],
                "refAccIntronStart": refSpliceInfo["intronStart"],
                "refDeNovoAccMES": deNovoAccInfo["refMaxEntScanScore"],
                "refDeNovoAccZ": deNovoAccInfo["refZScore"],
                "altDeNovoAccMES": deNovoAccInfo["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccInfo["altZScore"],
                "refDeNovoAccSeq": deNovoAccInfo["refSeq"],
                "altDeNovoAccSeq": deNovoAccInfo["altSeq"],
                "deNovoAccVarStart": deNovoAccInfo["varStart"],
                "deNovoAccVarLength": deNovoAccInfo["varLength"],
                "deNovoAccExonStart": deNovoAccInfo["exonStart"],
                "deNovoAccIntronStart": deNovoAccInfo["intronStart"],
                "closestAccRefMES": deNovoAccInfo["closestRefMaxEntScanScore"],
                "closestAccRefZ": deNovoAccInfo["closestRefZScore"],
                "closestAccRefSeq": deNovoAccInfo["closestRefSeq"],
                "closestAccAltMES": deNovoAccInfo["closestAltMaxEntScanScore"],
                "closestAccAltZ": deNovoAccInfo["closestAltZScore"],
                "closestAccAltSeq": deNovoAccInfo["closestAltSeq"],
                "closestAccExonStart": deNovoAccInfo["closestExonStart"],
                "closestAccIntronStart": deNovoAccInfo["closestIntronStart"],
                "deNovoAccAltGreaterRefFlag": deNovoAccInfo["altGreaterRefFlag"],
                "deNovoAccAltGreaterClosestRefFlag": deNovoAccInfo["altGreaterClosestRefFlag"],
                "deNovoAccAltGreaterClosestAltFlag": deNovoAccInfo["altGreaterClosestAltFlag"],
                "deNovoAccFrameshiftFlag": deNovoAccInfo["frameshiftFlag"],
                "spliceSite": refSpliceInfo["spliceSite"],
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag}
    
def getPriorProbProteinSNS(variant, variantData):
    '''
    Given a variant and a list of dictionaries containing variant data,
    Returns a dictionary containing:
      the variant's protein prior probability and enigma class for that prior
    '''
    proteinPrior = "-"
    enigmaClass = "-"
    if getVarType(variant) == "substitution":
        # TO DO possibly change HGVS_cDNA to pyhgvs_cDNA
        # [12:] parses out the NM_ accession so varHGVS is in format c.65C>T
        #varHGVS = variant["pyhgvs_cDNA"][12:]
        varHGVS = variant["HGVS_cDNA"]
        varGene = variant["Gene_Symbol"]

        for var in variantData:
            if var['gene'] == varGene and var['nthgvs'] == varHGVS:
                proteinPrior = float(var["protein_prior"])
                enigmaClass = getEnigmaClass(proteinPrior)
            
        return {"priorProb": proteinPrior,
                "enigmaClass": enigmaClass}

def getPriorProbInGreyZoneSNS(variant, boundaries, variantData):
    '''
    Given a variant and a list of dicitionaries with variant data,
    Returns applicable prior and enigma class based on protein priors for that variant
    Dictionary also contains other values that are either "N/A", "-", or 0 because they are not relevant
    '''
    if getVarType(variant) == "substitution" and getVarLocation(variant, boundaries) == "grey_zone_variant":
        proteinData = getPriorProbProteinSNS(variant, variantData)
        proteinPrior = proteinData["priorProb"]
        if proteinPrior == 0.99:
            proteinPrior = 0.5

        return {"applicablePrior": proteinPrior,
                "applicableEnigmaClass": getEnigmaClass(proteinPrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
                "refRefDonorMES": "-",
                "refRefDonorZ": "-",
                "altRefDonorMES": "-",
                "altRefDonorZ": "-",
                "refRefDonorSeq": "-",
                "altRefDonorSeq": "-",
                "refDonorVarStart": "-",
                "refDonorVarLength": "-",
                "refDonorExonStart": "-",
                "refDonorIntronStart": "-",
                "refDeNovoDonorMES": "-",
                "refDeNovoDonorZ": "-",
                "altDeNovoDonorMES": "-",
                "altDeNovoDonorZ": "-",
                "refDeNovoDonorSeq": "-",
                "altDeNovoDonorSeq": "-",
                "deNovoDonorVarStart": "-",
                "deNovoDonorVarLength": "-",
                "deNovoDonorExonStart": "-",
                "deNovoDonorIntronStart": "-",
                "closestDonorRefMES": "-",
                "closestDonorRefZ": "-",
                "closestDonorRefSeq": "-",
                "closestDonorAltMES": "-",
                "closestDonorAltZ": "-",
                "closestDonorAltSeq": "-",
                "closestDonorExonStart": "-",
                "closestDonorIntronStart": "-",
                "deNovoDonorAltGreaterRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestAltFlag": "N/A",
                "deNovoDonorFrameshiftFlag": "N/A",
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
                "refRefAccMES": "-",
                "refRefAccZ": "-",
                "altRefAccMES": "-",
                "altRefAccZ": "-",
                "refRefAccSeq": "-",
                "altRefAccSeq": "-",
                "refAccVarStart": "-",
                "refAccVarLength": "-",
                "refAccExonStart": "-",
                "refAccIntronStart": "-",
                "refDeNovoAccMES": "-",
                "refDeNovoAccZ": "-",
                "altDeNovoAccMES": "-",
                "altDeNovoAccZ": "-",
                "refDeNovoAccSeq": "-",
                "altDeNovoAccSeq": "-",
                "deNovoAccVarStart": "-",
                "deNovoAccVarLength": "-",
                "deNovoAccExonStart": "-",
                "deNovoAccIntronStart": "-",
                "closestAccRefMES": "-",
                "closestAccRefZ": "-",
                "closestAccRefSeq": "-",
                "closestAccAltMES": "-",
                "closestAccAltZ": "-",
                "closestAccAltSeq": "-",
                "closestAccExonStart": "-",
                "closestAccIntronStart": "-",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}
    
def getPriorProbInExonSNS(variant, boundaries, variantData):
    '''
    Given a variant, boundaries (either "enigma" or "priors") and a list of dictionaries containing variant data:
      1. Checks that variant is in an exon or clinically important domains and NOT in a splice site
      2. Checks that variant is a SNS variant
      3. Gets protein prior from variantData
      4. Determines if variant is a nonsense variant, if yes determines if splice rescue occurs
      5. If not a nonsense variant, calculates de novo donor prior and de novo acceptor prior if applicable
         Gets applicable prior if variant has a de novo donor prior
    Returns a dictionary containing all values, dictionary entry is "-" if not relevant to variant
    Values in dictionary include:
        applicable prior, highest prior if variant has multiple priors
        applicable class, highest predicted qualitative enigma class
        splice site = 0 because these variants are not in reference splice sites
        de novo donor and acceptor priors, prior prob for de novo splice sequence
        de novo acc flags (altGreaterRef, altGreaterClosestRef, altGreaterClosestAlt) = 1 if variant is possible de novo acceptor variant
        de novo donor flags (altGreaterRef, altGreaterClosestRef, altGreaterClosestAlt) = 1 if variant is possible de novo donor variant
        de novo acc/donor frameshift flags = 1 if de novo donor or acceptor causes a frameshfit
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshiftFlag = 1 if nonsense variant causes a frameshift mutation also, 0 if not, N/A if not nonsense variant
        inExonicPortionFlag = 1 if variant IS in exonic portion of highest scoring window and variant is a nonsense variant, 
           0 if NOT in exonic portion, N/A if not nonsense variant
        CIDomainRegionFlag = 1 if truncating region (between new donor and next reference splice acceptor) 
           includes a clinically import domain for nonsense variant, 0 if not, N/A if not nonsense variant
        isDivisibleFlag = 1 if distance between de novo donor and wild-type donor splice site is NOT divisible by 3 for nonsense variant, 
           0 if it IS divisible, N/A if not a nonsense variant
    '''
    varLoc = getVarLocation(variant, boundaries)
    if (varLoc == "exon_variant" or "CI_domain_variant") and getVarType(variant) == "substitution":
        proteinData = getPriorProbProteinSNS(variant, variantData)
        spliceRescue = "N/A"
        spliceFlag = 0
        frameshiftFlag = "N/A"
        inExonicPortionFlag = "N/A"
        CIDomainInRegionFlag = "N/A"
        isDivisibleFlag = "N/A"
        lowMESFlag = "N/A"
        deNovoDonorData = getPriorProbDeNovoDonorSNS(variant, STD_EXONIC_PORTION, deNovoDonorInRefAcc=False)
        if varInSpliceRegion(variant, donor=False, deNovo=True):
            deNovoAccData = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        else:
            deNovoAccData = {"priorProb": "N/A",
                             "refMaxEntScanScore": "N/A",
                             "altMaxEntScanScore": "N/A",
                             "refZScore": "N/A",
                             "altZScore": "N/A",
                             "refSeq": "N/A",
                             "altSeq": "N/A",
                             "varStart": "N/A",
                             "varLength": "N/A",
                             "exonStart": "N/A",
                             "intronStart": "N/A",
                             "closestRefMaxEntScanScore": "N/A",
                             "closestRefZScore": "N/A",
                             "closestRefSeq": "N/A",
                             "closestAltMaxEntScanScore": "N/A",
                             "closestAltZScore": "N/A",
                             "closestAltSeq": "N/A",
                             "closestExonStart": "N/A",
                             "closestIntronStart": "N/A",
                             "altGreaterRefFlag": "N/A",
                             "altGreaterClosestRefFlag": "N/A",
                             "altGreaterClosestAltFlag": "N/A",
                             "frameshiftFlag": "N/A"}
        varCons = getVarConsequences(variant)
        if varCons == "stop_gained":
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=False)
            applicablePrior = nonsenseData["priorProb"]
            applicableClass = nonsenseData["enigmaClass"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshiftFlag = nonsenseData["frameshiftFlag"]
            inExonicPortionFlag = nonsenseData["inExonicPortionFlag"]
            CIDomainInRegionFlag = nonsenseData["CIDomainInRegionFlag"]
            isDivisibleFlag = nonsenseData["isDivisibleFlag"]
            lowMESFlag = nonsenseData["lowMESFlag"]
        else:
            applicablePrior = proteinData["priorProb"]
            applicableClass = proteinData["enigmaClass"]
            if deNovoDonorData["priorProb"] != "N/A":
                applicablePrior = max(proteinData["priorProb"], deNovoDonorData["priorProb"])
                applicableClass = getEnigmaClass(applicablePrior)
                
        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": applicableClass,
                "proteinPrior": proteinData["priorProb"],
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": deNovoDonorData["priorProb"],
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorData["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorData["refZScore"],
                "altDeNovoDonorMES": deNovoDonorData["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorData["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorData["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorData["altSeq"],
                "deNovoDonorVarStart": deNovoDonorData["varStart"],
                "deNovoDonorVarLength": deNovoDonorData["varLength"],
                "deNovoDonorExonStart": deNovoDonorData["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorData["intronStart"],
                "closestDonorRefMES": deNovoDonorData["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorData["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorData["closestRefSeq"],
                "closestDonorAltMES": deNovoDonorData["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoDonorData["closestAltZScore"],
                "closestDonorAltSeq": deNovoDonorData["closestAltSeq"],
                "closestDonorExonStart": deNovoDonorData["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorData["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorData["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorData["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorData["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorData["frameshiftFlag"],
                "deNovoAccPrior": deNovoAccData["priorProb"],
                "refAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": deNovoAccData["refMaxEntScanScore"], 
                "refDeNovoAccZ": deNovoAccData["refZScore"],
                "altDeNovoAccMES": deNovoAccData["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccData["altZScore"],
                "refDeNovoAccSeq": deNovoAccData["refSeq"],
                "altDeNovoAccSeq": deNovoAccData["altSeq"],
                "deNovoAccVarStart": deNovoAccData["varStart"],
                "deNovoAccVarLength": deNovoAccData["varLength"],
                "deNovoAccExonStart": deNovoAccData["exonStart"],
                "deNovoAccIntronStart": deNovoAccData["intronStart"],
                "closestAccRefMES": deNovoAccData["closestRefMaxEntScanScore"],
                "closestAccRefZ": deNovoAccData["closestRefZScore"],
                "closestAccRefSeq": deNovoAccData["closestRefSeq"],
                "closestAccAltMES": deNovoAccData["closestAltMaxEntScanScore"],
                "closestAccAltZ": deNovoAccData["closestAltZScore"],
                "closestAccAltSeq": deNovoAccData["closestAltSeq"],
                "closestAccExonStart": deNovoAccData["closestExonStart"],
                "closestAccIntronStart": deNovoAccData["closestIntronStart"],
                "deNovoAccAltGreaterRefFlag": deNovoAccData["altGreaterRefFlag"],
                "deNovoAccAltGreaterClosestRefFlag": deNovoAccData["altGreaterClosestRefFlag"],
                "deNovoAccAltGreaterClosestAltFlag": deNovoAccData["altGreaterClosestAltFlag"],
                "deNovoAccFrameshiftFlag": deNovoAccData["frameshiftFlag"],
                "spliceSite": 0,
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshiftFlag": frameshiftFlag,
                "inExonicPortionFlag": inExonicPortionFlag,
                "CIDomainInRegionFlag": CIDomainInRegionFlag,
                "isDivisibleFlag": isDivisibleFlag,
                "lowMESFlag": lowMESFlag}

def getPriorProbOutsideTranscriptBoundsSNS(variant, boundaries):
    '''
    Given a variant and boundaries (either "enigma" or "priors"),
    Checks that variant is outside transcript boundaries
    Returns prior prob and predicted qualitative enigma class
    Dictionary also contains other values that are either "N/A", "-", or 0 because they are not relevant
    '''
    varLoc = getVarLocation(variant, boundaries)
    varType = getVarType(variant)
    if varLoc == "outside_transcript_boundaries_variant" and varType == "substitution":
        priorProb = 0.02
        return {"applicablePrior": priorProb,
                "applicableEnigmaClass": getEnigmaClass(priorProb),
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
                "refRefDonorMES": "-",
                "refRefDonorZ": "-",
                "altRefDonorMES": "-",
                "altRefDonorZ": "-",
                "refRefDonorSeq": "-",
                "altRefDonorSeq": "-",
                "refDonorVarStart": "-",
                "refDonorVarLength": "-",
                "refDonorExonStart": "-",
                "refDonorIntronStart": "-",
                "refDeNovoDonorMES": "-",
                "refDeNovoDonorZ": "-",
                "altDeNovoDonorMES": "-",
                "altDeNovoDonorZ": "-",
                "refDeNovoDonorSeq": "-",
                "altDeNovoDonorSeq": "-",
                "deNovoDonorVarStart": "-",
                "deNovoDonorVarLength": "-",
                "deNovoDonorExonStart": "-",
                "deNovoDonorIntronStart": "-",
                "closestDonorRefMES": "-",
                "closestDonorRefZ": "-",
                "closestDonorRefSeq": "-",
                "closestDonorAltMES": "-",
                "closestDonorAltZ": "-",
                "closestDonorAltSeq": "-",
                "closestDonorExonStart": "-",
                "closestDonorIntronStart": "-",
                "deNovoDonorAltGreaterRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestRefFlag": "N/A",
                "deNovoDonorAltGreaterClosestAltFlag": "N/A",
                "deNovoDonorFrameshiftFlag": "N/A",
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
                "refRefAccMES": "-",
                "refRefAccZ": "-",
                "altRefAccMES": "-",
                "altRefAccZ": "-",
                "refRefAccSeq": "-",
                "altRefAccSeq": "-",
                "refAccVarStart": "-",
                "refAccVarLength": "-",
                "refAccExonStart": "-",
                "refAccIntronStart": "-",
                "refDeNovoAccMES": "-",
                "refDeNovoAccZ": "-",
                "altDeNovoAccMES": "-",
                "altDeNovoAccZ": "-",
                "refDeNovoAccSeq": "-",
                "altDeNovoAccSeq": "-",
                "deNovoAccVarStart": "-",
                "deNovoAccVarLength": "-",
                "deNovoAccExonStart": "-",
                "deNovoAccIntronStart": "-",
                "closestAccRefMES": "-",
                "closestAccRefZ": "-",
                "closestAccRefSeq": "-",
                "closestAccAltMES": "-",
                "closestAccAltZ": "-",
                "closestAccAltSeq": "-",
                "closestAccExonStart": "-",
                "closestAccIntronStart": "-",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": 0,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}

def getPriorProbIntronicDeNovoDonorSNS(variant):
    '''
    Given a variant,
      1. Checks that variant is NOT in exon or reference donor/acceptor site
      2. Checks that variant is a substitution variant
    Determines if alt MES score is greater than ref MES score for highest scoring sliding window
      If altMES > refMES, then deNovoDonorAltGreaterRefFlag = 1 (0 otherwise)
    Determines closest ref donor scores
      If altMES > closestRefDonorMES, then deNovoDonorAltGreaterClosestRefFlag = 1 (0 otherwise)
    If either flag is equal to 1, flag variant for further analysis (spliceFlag = 1), spliceFlag = 0 otherwise
    altGreaterClosestAltFlag is always equals N/A because this function is not used for any variants in ref splice sites
    Returns dictionary containing prior prob, enigma class, de novo donor scores, and splice flag
    '''
    inExon = varInExon(variant)
    inRefDonor = varInSpliceRegion(variant, donor=True, deNovo=False)
    inRefAcc = varInSpliceRegion(variant, donor=False, deNovo=False)
    if inExon == False and inRefDonor == False and inRefAcc == False:
        if getVarType(variant) == "substitution":
            deNovoDonorInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                                    donor=True, deNovo=False, deNovoDonorInRefAcc=False)
            refMES = deNovoDonorInfo["refMaxEntScanScore"]
            altMES = deNovoDonorInfo["altMaxEntScanScore"]
            deNovoOffset = 0
            closestDonorInfo = getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False, deNovoDonorInRefAcc=False, testMode=False)
            closestMES = closestDonorInfo["maxEntScanScore"]
            spliceFlag = 0
            altGreaterRefFlag = 0
            altGreaterClosestRefFlag = 0
            if altMES > refMES:
                spliceFlag = 1
                altGreaterRefFlag = 1
            if altMES > closestMES:
                spliceFlag = 1
                altGreaterClosestRefFlag = 1

            frameshiftStatus = getDeNovoSpliceFrameshiftStatus(variant, donor=True, deNovoDonorInRefAcc=False)
            frameshiftFlag = 0
            if frameshiftStatus == True:
                frameshiftFlag = 1

            return {"priorProb": "N/A",
                    "enigmaClass": "N/A",
                    "refMaxEntScanScore": refMES,
                    "refZScore": deNovoDonorInfo["refZScore"],
                    "altMaxEntScanScore": altMES,
                    "altZScore": deNovoDonorInfo["altZScore"],
                    "refSeq": deNovoDonorInfo["refSeq"],
                    "altSeq": deNovoDonorInfo["altSeq"],
                    "varStart": deNovoDonorInfo["varStart"],
                    "varLength": deNovoDonorInfo["varLength"],
                    "exonStart": 0,
                    "intronStart": STD_EXONIC_PORTION,
                    "closestRefSeq": closestDonorInfo["sequence"],
                    "closestRefMaxEntScanScore": closestMES,
                    "closestRefZScore": closestDonorInfo["zScore"],
                    "closestAltSeq": "N/A",
                    "closestAltMaxEntScanScore": "N/A",
                    "closestAltZScore": "N/A",
                    "closestExonStart": closestDonorInfo["exonStart"],
                    "closestIntronStart": closestDonorInfo["intronStart"],
                    "altGreaterRefFlag": altGreaterRefFlag,
                    "altGreaterClosestRefFlag": altGreaterClosestRefFlag,
                    "altGreaterClosestAltFlag": "N/A",
                    "frameshiftFlag": frameshiftFlag,
                    "spliceFlag": spliceFlag}
    
def getPriorProbInIntronSNS(variant, boundaries):
    '''
    Given a variant and boundaries (either "priors or "enigma"),
    Checks that variant is located in an intron and is a substitution variant
    Determines if variant creates a de novo donor site in the intron
    Returns a dictionary containing applicable prior and predicted qualitative enigma class
    Dictionary also contains de novo donor ref and alt scores
    AND a spliceFlag which is equal to 1 if variant creates a better de novo splice site than ref sequence
    Rest of values in dictionary are equal to 0, "-", or N/A because they are not relevant to variant
    '''
    varLoc = getVarLocation(variant, boundaries)
    varType = getVarType(variant)
    if varLoc == "intron_variant" and varType == "substitution":
        deNovoDonorData = getPriorProbIntronicDeNovoDonorSNS(variant)
        if deNovoDonorData["spliceFlag"] == 0:
            priorProb = 0.02
            enigmaClass = getEnigmaClass(priorProb)
        else:
            priorProb = deNovoDonorData["priorProb"]
            enigmaClass = deNovoDonorData["enigmaClass"]
            
        return {"applicablePrior": priorProb,
                "applicableEnigmaClass": enigmaClass,
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": "N/A",
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorData["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorData["refZScore"],
                "altDeNovoDonorMES": deNovoDonorData["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorData["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorData["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorData["altSeq"],
                "deNovoDonorVarStart": deNovoDonorData["varStart"],
                "deNovoDonorVarLength": deNovoDonorData["varLength"],
                "deNovoDonorExonStart": deNovoDonorData["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorData["intronStart"],
                "closestDonorRefMES": deNovoDonorData["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorData["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorData["closestRefSeq"],
                "closestDonorAltMES": "N/A",
                "closestDonorAltZ": "N/A",
                "closestDonorAltSeq": "N/A",
                "closestDonorExonStart": deNovoDonorData["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorData["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorData["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorData["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorData["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorData["frameshiftFlag"],
                "deNovoAccPrior": "N/A",
                "refAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",
                "refDeNovoAccMES": "N/A",
                "refDeNovoAccZ": "N/A",
                "altDeNovoAccMES": "N/A",
                "altDeNovoAccZ": "N/A",
                "refDeNovoAccSeq": "N/A",
                "altDeNovoAccSeq": "N/A",
                "deNovoAccVarStart": "N/A",
                "deNovoAccVarLength": "N/A",
                "deNovoAccExonStart": "N/A",
                "deNovoAccIntronStart": "N/A",
                "closestAccRefMES": "N/A",
                "closestAccRefZ": "N/A",
                "closestAccRefSeq": "N/A",
                "closestAccAltMES": "N/A",
                "closestAccAltZ": "N/A",
                "closestAccAltSeq": "N/A",
                "closestAccExonStart": "N/A",
                "closestAccIntronStart": "N/A",
                "deNovoAccAltGreaterRefFlag": "N/A",
                "deNovoAccAltGreaterClosestRefFlag": "N/A",
                "deNovoAccAltGreaterClosestAltFlag": "N/A",
                "deNovoAccFrameshiftFlag": "N/A",
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": deNovoDonorData["spliceFlag"],
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}

def getPriorProbInUTRSNS(variant, boundaries):
    '''
    Given a variant and boundaries (either "priors" or "enigma"),
    Checks that variant is a SNS variant in a UTR
    Determines prior prob based on location (5'/3' UTR and intron/exon)
    Returns a dictionary containing applicable prior and predicted qualitative enigma class
    Dictionary also contains de novo donor/acceptor ref and alt scores if applicable
    AND a spliceFlag which is equal to 1 if variant creates a better de novo splice site than ref sequence
    Rest of values in dictionary are equal to 0, "-", or N/A because they are not relevant to variant
    '''
    # TO DO still need to account for creation of alternate ATG codons in 5' UTRs
    varLoc = getVarLocation(variant, boundaries)
    varType = getVarType(variant)
    if varLoc == "UTR_variant" and varType == "substitution":
        deNovoAccData = {"priorProb": "N/A",
                         "refMaxEntScanScore": "N/A",
                         "altMaxEntScanScore": "N/A",
                         "refZScore": "N/A",
                         "altZScore": "N/A",
                         "refSeq": "N/A",
                         "altSeq": "N/A",
                         "varStart": "N/A",
                         "varLength": "N/A",
                         "exonStart": "N/A",
                         "intronStart": "N/A",
                         "closestRefMaxEntScanScore": "N/A",
                         "closestRefZScore": "N/A",
                         "closestRefSeq": "N/A",
                         "closestAltMaxEntScanScore": "N/A",
                         "closestAltZScore": "N/A",
                         "closestAltSeq": "N/A",
                         "closestExonStart": "N/A",
                         "closestIntronStart": "N/A",
                         "altGreaterRefFlag": "N/A",
                         "altGreaterClosestRefFlag": "N/A",
                         "altGreaterClosestAltFlag": "N/A",
                         "frameshiftFlag": "N/A"}
        deNovoDonorData = {"priorProb": "N/A",
                           "refMaxEntScanScore": "N/A",
                           "altMaxEntScanScore": "N/A",
                           "refZScore": "N/A",
                           "altZScore": "N/A",
                           "refSeq": "N/A",
                           "altSeq": "N/A",
                           "varStart": "N/A",
                           "varLength": "N/A",
                           "exonStart": "N/A",
                           "intronStart": "N/A",
                           "closestRefMaxEntScanScore": "N/A",
                           "closestRefZScore": "N/A",
                           "closestRefSeq": "N/A",
                           "closestAltMaxEntScanScore": "N/A",
                           "closestAltZScore": "N/A",
                           "closestAltSeq": "N/A",
                           "closestExonStart": "N/A",
                           "closestIntronStart": "N/A",
                           "altGreaterRefFlag": "N/A",
                           "altGreaterClosestRefFlag": "N/A",
                           "altGreaterClosestAltFlag": "N/A",
                           "frameshiftFlag": "N/A"}
        varCons = getVarConsequences(variant)
        if varCons == "3_prime_UTR_variant":
            applicablePrior = 0.02
            applicableClass = getEnigmaClass(applicablePrior)
            spliceFlag = 0
        elif varCons == "5_prime_UTR_variant":
            if varInExon(variant) == True:
                deNovoDonorData = getPriorProbDeNovoDonorSNS(variant, STD_EXONIC_PORTION, deNovoDonorInRefAcc=False)
                applicablePrior = deNovoDonorData["priorProb"]
                applicableClass = deNovoDonorData["enigmaClass"]
                spliceFlag = 0
                if varInSpliceRegion(variant, donor=False, deNovo=True) == True:
                    deNovoAccData = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        else:
            # to account for variants in 5' UTR that are classified as other variant types
            deNovoDonorData = getPriorProbIntronicDeNovoDonorSNS(variant)
            spliceFlag = deNovoDonorData["spliceFlag"]
            if spliceFlag == 1:
                applicablePrior = deNovoDonorData["priorProb"]
                applicableClass = deNovoDonorData["enigmaClass"]
            else:
                applicablePrior = 0.02
                applicableClass = getEnigmaClass(applicablePrior)

        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": applicableClass,
                "proteinPrior": "N/A",
                "refDonorPrior": "N/A",
                "deNovoDonorPrior": deNovoDonorData["priorProb"],
                "refRefDonorMES": "N/A",
                "refRefDonorZ": "N/A",
                "altRefDonorMES": "N/A",
                "altRefDonorZ": "N/A",
                "refRefDonorSeq": "N/A",
                "altRefDonorSeq": "N/A",
                "refDonorVarStart": "N/A",
                "refDonorVarLength": "N/A",
                "refDonorExonStart": "N/A",
                "refDonorIntronStart": "N/A",
                "refDeNovoDonorMES": deNovoDonorData["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorData["refZScore"],
                "altDeNovoDonorMES": deNovoDonorData["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorData["altZScore"],
                "refDeNovoDonorSeq": deNovoDonorData["refSeq"],
                "altDeNovoDonorSeq": deNovoDonorData["altSeq"],
                "deNovoDonorVarStart": deNovoDonorData["varStart"],
                "deNovoDonorVarLength": deNovoDonorData["varLength"],
                "deNovoDonorExonStart": deNovoDonorData["exonStart"],
                "deNovoDonorIntronStart": deNovoDonorData["intronStart"],
                "closestDonorRefMES": deNovoDonorData["closestRefMaxEntScanScore"],
                "closestDonorRefZ": deNovoDonorData["closestRefZScore"],
                "closestDonorRefSeq": deNovoDonorData["closestRefSeq"],
                "closestDonorAltMES": deNovoDonorData["closestAltMaxEntScanScore"],
                "closestDonorAltZ": deNovoDonorData["closestAltZScore"],
                "closestDonorAltSeq": deNovoDonorData["closestAltSeq"],
                "closestDonorExonStart": deNovoDonorData["closestExonStart"],
                "closestDonorIntronStart": deNovoDonorData["closestIntronStart"],
                "deNovoDonorAltGreaterRefFlag": deNovoDonorData["altGreaterRefFlag"],
                "deNovoDonorAltGreaterClosestRefFlag": deNovoDonorData["altGreaterClosestRefFlag"],
                "deNovoDonorAltGreaterClosestAltFlag": deNovoDonorData["altGreaterClosestAltFlag"],
                "deNovoDonorFrameshiftFlag": deNovoDonorData["frameshiftFlag"],
                "deNovoAccPrior": deNovoAccData["priorProb"],
                "refAccPrior": "N/A",
                "refRefAccMES": "N/A",
                "refRefAccZ": "N/A",
                "altRefAccMES": "N/A",
                "altRefAccZ": "N/A",
                "refRefAccSeq": "N/A",
                "altRefAccSeq": "N/A",
                "refAccVarStart": "N/A",
                "refAccVarLength": "N/A",
                "refAccExonStart": "N/A",
                "refAccIntronStart": "N/A",                
                "refDeNovoAccMES": deNovoAccData["refMaxEntScanScore"], 
                "refDeNovoAccZ": deNovoAccData["refZScore"],
                "altDeNovoAccMES": deNovoAccData["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccData["altZScore"],
                "refDeNovoAccSeq": deNovoAccData["refSeq"],
                "altDeNovoAccSeq": deNovoAccData["altSeq"],
                "deNovoAccVarStart": deNovoAccData["varStart"],
                "deNovoAccVarLength": deNovoAccData["varLength"],
                "deNovoAccExonStart": deNovoAccData["exonStart"],
                "deNovoAccIntronStart": deNovoAccData["exonStart"],
                "closestAccRefMES": deNovoAccData["closestRefMaxEntScanScore"],
                "closestAccRefZ": deNovoAccData["closestRefZScore"],
                "closestAccRefSeq": deNovoAccData["closestRefSeq"],
                "closestAccAltMES": deNovoAccData["closestAltMaxEntScanScore"],
                "closestAccAltZ": deNovoAccData["closestAltZScore"],
                "closestAccAltSeq": deNovoAccData["closestAltSeq"],
                "closestAccExonStart": deNovoAccData["closestExonStart"],
                "closestAccIntronStart": deNovoAccData["closestExonStart"],
                "deNovoAccAltGreaterRefFlag": deNovoAccData["altGreaterRefFlag"],
                "deNovoAccAltGreaterClosestRefFlag": deNovoAccData["altGreaterClosestRefFlag"],
                "deNovoAccAltGreaterClosestAltFlag": deNovoAccData["altGreaterClosestAltFlag"],
                "deNovoAccFrameshiftFlag": deNovoAccData["frameshiftFlag"],
                "spliceSite": 0,
                "spliceRescue": "N/A",
                "spliceFlag": spliceFlag,
                "frameshiftFlag": "N/A",
                "inExonicPortionFlag": "N/A",
                "CIDomainInRegionFlag": "N/A",
                "isDivisibleFlag": "N/A",
                "lowMESFlag": "N/A"}

def getVarData(variant, boundaries, variantData):
    '''
    Given variant, boundaries (either "priors" or "enigma') and list of dictionaries with variant data
    Checks that variant is a single nucleotide substitution
    Determines prior prob dictionary based on variant location
    Return dictionary containing all values for all new prior prob fields
    '''
    varLoc = getVarLocation(variant, boundaries)
    varType = getVarType(variant)
    blankDict = {"applicablePrior": "-",
                 "applicableEnigmaClass": "-",
                 "proteinPrior": "-",
                 "refDonorPrior": "-",
                 "deNovoDonorPrior": "-",
                 "refRefDonorMES": "-",
                 "refRefDonorZ": "-",
                 "altRefDonorMES": "-",
                 "altRefDonorZ": "-",
                 "refRefDonorSeq": "-",
                 "altRefDonorSeq": "-",
                 "refDonorVarStart": "-",
                 "refDonorVarLength": "-",
                 "refDonorExonStart": "-",
                 "refDonorIntronStart": "-",                
                 "refDeNovoDonorMES": "-",
                 "refDeNovoDonorZ": "-",
                 "altDeNovoDonorMES": "-",
                 "altDeNovoDonorZ": "-",
                 "refDeNovoDonorSeq": "-",
                 "altDeNovoDonorSeq": "-",
                 "deNovoDonorVarStart": "-",
                 "deNovoDonorVarLength": "-",
                 "deNovoDonorExonStart": "-",
                 "deNovoDonorIntronStart": "-",
                 "closestDonorRefMES": "-",
                 "closestDonorRefZ": "-",
                 "closestDonorRefSeq": "-",
                 "closestDonorAltMES": "-",
                 "closestDonorAltZ": "-",
                 "closestDonorAltSeq": "-",
                 "closestDonorExonStart": "-",
                 "closestDonorIntronStart": "-",
                 "deNovoDonorAltGreaterRefFlag": "-",
                 "deNovoDonorAltGreaterClosestRefFlag": "-",
                 "deNovoDonorAltGreaterClosestAltFlag": "-",
                 "deNovoDonorFrameshiftFlag": "-",
                 "refAccPrior": "-",
                 "deNovoAccPrior": "-",
                 "refRefAccMES": "-",
                 "refRefAccZ": "-",
                 "altRefAccMES": "-",
                 "altRefAccZ": "-",
                 "refRefAccSeq": "-",
                 "altRefAccSeq": "-",
                 "refAccVarStart": "-",
                 "refAccVarLength": "-",
                 "refAccExonStart": "-",
                 "refAccIntronStart": "-",                 
                 "refDeNovoAccMES": "-",
                 "refDeNovoAccZ": "-",
                 "altDeNovoAccMES": "-",
                 "altDeNovoAccZ": "-",
                 "refDeNovoAccSeq": "-",
                 "altDeNovoAccSeq": "-",
                 "deNovoAccVarStart": "-",
                 "deNovoAccVarLength": "-",
                 "deNovoAccExonStart": "-",
                 "deNovoAccIntronStart": "-",
                 "closestAccRefMES": "-",
                 "closestAccRefZ": "-",
                 "closestAccRefSeq": "-",
                 "closestAccAltMES": "-",
                 "closestAccAltZ": "-",
                 "closestAccAltSeq": "-",
                 "closestAccExonStart": "-",
                 "closestAccIntronStart": "-",
                 "deNovoAccAltGreaterRefFlag": "-",
                 "deNovoAccAltGreaterClosestRefFlag": "-",
                 "deNovoAccAltGreaterClosestAltFlag": "-",
                 "deNovoAccFrameshiftFlag": "-",
                 "spliceSite": "-",
                 "spliceRescue": "-",
                 "spliceFlag": "-",
                 "frameshiftFlag": "-",
                 "inExonicPortionFlag": "-",
                 "CIDomainInRegionFlag": "-",
                 "isDivisibleFlag": "-",
                 "lowMESFlag": "-"}
    if varType == "substitution":
        # functions only work for variants with cannonical nucleotides (ACTG)
        if variant["Ref"] in ["A", "C","G", "T"] and variant["Alt"] in ["A", "C", "G", "T"]:
            if varLoc == "outside_transcript_boundaries_variant":
                varData = getPriorProbOutsideTranscriptBoundsSNS(variant, boundaries)
            elif varLoc == "CI_splice_donor_variant" or varLoc == "splice_donor_variant":
                varData =  getPriorProbSpliceDonorSNS(variant, boundaries, variantData)
            elif varLoc == "CI_splice_acceptor_variant" or varLoc == "splice_acceptor_variant":
                varData =  getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData)
            elif varLoc == "CI_domain_variant" or varLoc == "exon_variant":
                varData = getPriorProbInExonSNS(variant, boundaries, variantData)
            elif varLoc == "grey_zone_variant":
                varData = getPriorProbInGreyZoneSNS(variant, boundaries, variantData)
            elif varLoc == "after_grey_zone_variant":
                varData =  getPriorProbAfterGreyZoneSNS(variant, boundaries)
            elif varLoc == "UTR_variant":
                varData =  getPriorProbInUTRSNS(variant, boundaries)
            elif varLoc == "intron_variant":
                varData =  getPriorProbInIntronSNS(variant, boundaries)
            else:
                varData = blankDict.copy()
        else:
            varData = blankDict.copy()
    else:
        # to account for any non SNS variants
        varData = blankDict.copy()
        
    varData["varType"] = varType
    varData["varLoc"] = varLoc
    return varData
            
def addVarDataToRow(varData, inputRow):
    '''
    Given data about a particular variant and a row from input file,
    Returns row with appended data
    '''
    for key in varData.keys():
        inputRow[key] = varData[key]
    return inputRow

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
               "varHGVScDNA": variant["HGVS_cDNA"]}

    return varDict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputFile", default="built.tsv", help="File with variant information")
    parser.add_argument('-o', "--outputFile", help="File where results will be output")
    parser.add_argument('-v', "--variantFile", help="File containing protein priors for variants")
    parser.add_argument('-b', "--boundaries", default="enigma",
                        help="Specifies which boundaries ('enigma' or 'priors') to use for clinically important domains")
    args = parser.parse_args()    

    inputData = csv.DictReader(open(args.inputFile, "r"), delimiter="\t")
    fieldnames = inputData.fieldnames
    newHeaders = ["varType", "varLoc", "applicablePrior", "applicableEnigmaClass", "proteinPrior", "refDonorPrior", "deNovoDonorPrior",
                  "refRefDonorMES", "refRefDonorZ", "altRefDonorMES", "altRefDonorZ", "refRefDonorSeq", "altRefDonorSeq", "refDonorVarStart",
                  "refDonorVarLength", "refDonorExonStart", "refDonorIntronStart", "refDeNovoDonorMES", "refDeNovoDonorZ", "altDeNovoDonorMES",
                  "altDeNovoDonorZ", "refDeNovoDonorSeq", "altDeNovoDonorSeq", "deNovoDonorVarStart", "deNovoDonorVarLength",
                  "deNovoDonorExonStart", "deNovoDonorIntronStart", "closestDonorRefMES", "closestDonorRefZ", "closestDonorRefSeq",
                  "closestDonorAltMES", "closestDonorAltZ", "closestDonorAltSeq", "closestDonorExonStart",
                  "closestDonorIntronStart", "deNovoDonorAltGreaterRefFlag", "deNovoDonorAltGreaterClosestRefFlag",
                  "deNovoDonorAltGreaterClosestAltFlag", "deNovoDonorFrameshiftFlag", "refAccPrior", "deNovoAccPrior", "refRefAccMES",
                  "refRefAccZ", "altRefAccMES", "altRefAccZ", "refRefAccSeq", "altRefAccSeq", "refAccVarStart", "refAccVarLength", "refAccExonStart",
                  "refAccIntronStart", "refDeNovoAccMES", "refDeNovoAccZ", "altDeNovoAccMES", "altDeNovoAccZ", "refDeNovoAccSeq", "altDeNovoAccSeq",
                  "deNovoAccVarStart", "deNovoAccVarLength", "deNovoAccExonStart", "deNovoAccIntronStart",
                  "closestAccRefMES", "closestAccRefZ", "closestAccRefSeq", "closestAccAltMES", "closestAccAltZ", "closestAccAltSeq",
                  "closestAccExonStart", "closestAccIntronStart", "deNovoAccAltGreaterRefFlag", "deNovoAccAltGreaterClosestRefFlag",
                  "deNovoAccAltGreaterClosestAltFlag", "deNovoAccFrameshiftFlag", "spliceSite", "spliceRescue", "spliceFlag", "frameshiftFlag",
                  "inExonicPortionFlag", "CIDomainInRegionFlag", "isDivisibleFlag", "lowMESFlag"]
    for header in newHeaders:
        fieldnames.append(header)
    outputData = csv.DictWriter(open(args.outputFile, "w"), delimiter="\t", fieldnames=fieldnames)
    outputData.writerow(dict((fn,fn) for fn in inputData.fieldnames))

    totalVariants = 0
    for variant in inputData:
        variantData = csv.DictReader(open(args.variantFile, "r"), delimiter="\t")    
        varData = getVarData(variant, args.boundaries, variantData)
        variant = addVarDataToRow(varData, variant)
        outputData.writerow(variant)
        totalVariants += 1
        print "variant", totalVariants, "complete"
    
if __name__ == "__main__":
    main()
