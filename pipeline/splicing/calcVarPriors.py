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
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pygr.seqdb import SequenceFileDB
from Bio.Seq import Seq
from calcMaxEntScanMeanStd import fetch_gene_coordinates, runMaxEntScan

'''
GENERAL NOTES ON REFSEQ NUMBERING AND SPLICING

MINUS VS PLUS STRAND NUMBERING
Genes are transcribed from 5' to 3'
For plus strand genes (e.g. BRCA2) this correlates to reading from left to right
    so genomic position increases from left to right
This means that for plus strand genes the exon start position is ALWAYS less than
    the exon end position and that native splice acceptors have a lower genomic position
    than native splice donors in the same exon
For minus strand genes (e.g. BRCA1) this correlates to reading from right to left
    so genomic position increases from right to left AND decreases from left to right
This means that for minus strand genes the exon start position is ALWAYS greater than
    the exon end position and that native splice acceptors have a higher genomic position
    than native splice donors in the same exon

REFSEQ NUMBERING
RefSeq numbering starts to the right of the first base
This effects the numbering for the lefthand side of genes when reading from left to right
    in a layout where genomic position increases when you read from left to right
For plus strand genes (e.g. BRCA2) this means that the RefSeq exon start genomic position
    (which is the leftmost position in a layout where genomic position increases from left to right)
    refers to the base immediately to the left of the first base of the exon
    this base is the last base of the intron and is referred to as the -1 position in transcript coordinates 
        (e.g. c.476-1)
    so the exon actually starts at varExonStart + 1 which is the genomic position of the first base in the exon
        because numbering increases from left to right so this is the base to the left of varExonStart
For minus strand genes (e.g. BRCA1) this means that the RefSeq exon end genomic position position
    (which is the leftmost position in a layout where genomic position increases from left to right
        AND is the rightmost position in a layout where genomic position decreases form left to right)
    refers to the base immediately to the right of the last base of the exon (when looking at a layout
        where numbering decreases from left to right) 
    this base in the first base of the intorn and is referred to as the +1 position in transcript coordinates
    (e.g. 4185+1)
    so the exon actually ends at varExonEnd + 1 which is the genomic position of the last base in the exon
        because numbering decreases from left to right this is the base that is to the left of varExonEnd

SPLICING AND SPLICE POSITIONS
Native splice acceptors are located at the 5' ends of exons (near exonStart)
Native splice donors are located at the 3' ends of exons (near exonEnd)
For splice donors:
    splicing occurs between the last base of the exon and the first base of the intron
    this means that the last base in the exon is considered the wild-type splice cut position
        (the function to get FASTA sequences is inclusive of the endpoints so this is the last base included)
    this also means that the first base in the intron is considered the wild-type splice donor position
        (per HCI PRIORS website) 
    For example in a minus strand gene (BRCA1):
        c.441 (g.43104122) is the last base in the exon and is the wild-type cut position 
        and c.441+1 (g.43104121) is the first base in the intron and is the wild-type donor position
        So the donor position is to the right of the cut position 
            (when looking at a layout where numbering decreases from left to right)
    For example in a plus strand gene (BRCA2):
        c.7617 (g.32356609) is the last base in the exon and is the wild-type cut posiiton
        and c.7617+1 (g.32356610) is the first base in the intron and is the wild-type donor position
        So the donor position is to the right of the cut position
            (when looking at a layout where numbering increases from left to right)
For splice acceptors:
    splicing occurs between the last base of the intron and the first base of the exon
    this means that the first base in the exon is considered the wild-type splice cut position
        (the function to get FASTA sequences is inclusive of the endpoints to this is the first base included)
    this also means that the last base in the intron is considered the wild-type splice acceptor position
        (which matches the guidelines used for splice donors on the HCI PRIORS website)
    For example in a minus strand gene (BRCA1):
        c.5194 (g.43057135) is the first base in the exon and is the wild-type cut posiiton
        and c.5194-1 (g.43057136) is the last base in the intron and is the wild-type acceptor position
        So the acceptor position is to the left of the cut position
            (when looking at a layout where numbering decreases from left to right)
    For example in a plus strand gene (BRCA2):
        c.317 (g.32325076) is the first base in the exon and is the wild-type cut position
        and c.317-1 (g.32325075) is the last base in the intron and is the wild-type acceptor position
        So the acceptor position is to the left of the cut position
            (when looking at a layout where numbering increases from left to right)
'''

'''
NOTES ON PRIOR PROBABILITY DICTIONARIES

These dictionaries are generated by functions that start with "getPriorPrior"
All values in these dictionaries are assigned either numerical values, "N/A", or "-"
"N/A" is included when that particular dictionary key is not relevant for a given variant
"-" is included when a particular dictionary key is not evaluated (the value is not calculated)
'''

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

# Probability constants from SVT and MP valid as of 5/24/18
LOW_PROBABILITY = 0.02 
LOW_SPLICING_PROBABILITY = 0.04
MODERATE_DENOVO_PROBABILITY = 0.3
MODERATE_SPLICING_PROBABILITY = 0.34
CAPPED_PROBABILITY = 0.5
HIGH_PROBABILITY = 0.64
HIGH_SPLICING_PROBABILITY = 0.97
PATHOGENIC_PROBABILITY = 0.99

# Splicing cutoff constants from SVT and MP valid as of 5/24/18
LOW_MES_CUTOFF = 6.2
HIGH_MES_CUTOFF = 8.5
MIN_REF_ALT_ZDIFF = 0.5
REF_DONOR_REFZ_CUTOFF = -1.5
REF_DONOR_HIGH_CUTOFF = 0.0
REF_DONOR_LOW_CUTOFF = -2.0
REF_ACC_REFZ_CUTOFF = -1.0
REF_ACC_HIGH_CUTOFF = 0.5
REF_ACC_LOW_CUTOFF = -1.5
DE_NOVO_DONOR_HIGH_CUTOFF = 0.0
DE_NOVO_DONOR_LOW_CUTOFF = -2.0

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
                if len(varAlt) == 1 and varAlt == varRef[0]:
                    return "deletion"
                else:
                    return "delins"
            elif len(varRef) < len(varAlt):
                if len(varRef) == 1 and varRef == varAlt[0]:
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
                if varGenPos <= exonStart and varGenPos > exonEnd:
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
                if varGenPos <= exonStart and varGenPos > exonEnd:
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

def getFastaSeq(chrom, rangeStart, rangeStop, plusStrandSeq=True):
    '''
    Given chromosome (in format 'chr13'), region genomic start position, and
    region genomic end position:
    Returns a string containing the sequence inclusive of rangeStart and rangeStop
    If plusStrandSeq=True, returns plus strand sequence
    If plusStrandSeq=False, returns minus strand sequence
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
    for base in sequence:
        assert base in ["A", "C", "G", "T", "a", "c", "g", "t"]
    if plusStrandSeq == True:
        return sequence.upper()
    else:
        return str(Seq(sequence.upper()).reverse_complement())

def getSeqLocDict(chrom, varStrand, rangeStart, rangeStop):
    '''
    Given chromosome, strand, region genomic start position, and region genomic end position
    returns a dictionary containing the genomic position as the key and reference allele as the value
    For minus strand gene (rangeStart > rangeStop), for plus strand gene (rangeStart < rangeStop)
    Always returns plus strand sequence 
    '''
    seqLocDict = {}
    if varStrand == "-":
        regionStart = int(rangeStop)
        regionEnd = int(rangeStart)
    else:
        regionStart = int(rangeStart)
        regionEnd = int(rangeStop)
    sequence = getFastaSeq(chrom, regionStart, regionEnd, plusStrandSeq=True)
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

def getAltSeq(altSeqDict, varStrand):
    '''
    Given a dictionary containing an alternate sequence with bases and their locations
    and the strand that the alternate allele is on
    Returns a string of the sequence containing the alternate allele
    '''
    sequence = ""
    # to ensure that items in dictionary are sorted numerically
    for key, value in sorted(altSeqDict.items()):
        sequence += altSeqDict[key]
    if varStrand == "-":
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
        refSeq = getFastaSeq(varChrom, rangeStart, rangeStop, plusStrandSeq=False)
    else:
        refSeq = getFastaSeq(varChrom, rangeStart, rangeStop, plusStrandSeq=True)
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
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["donorStart"], refSpliceBounds["donorEnd"], plusStrandSeq=True)
            else:
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["donorStart"], refSpliceBounds["donorEnd"], plusStrandSeq=False)
        else:
            refSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=True)
            deNovoOffset = deNovoLength - exonicPortionSize
            # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            if getVarStrand(variant) == "+":
                # acceptorEnd - deNovoOffset because genomic position increases from left to right on plus strand, refSeq reduced to correct length
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] - deNovoOffset), plusStrandSeq=True)
            else:
                # acceptorEnd + deNovoOffset because genomic position decreases from left to right on minus strand, refSeq reduced to correct length
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] + deNovoOffset), plusStrandSeq=False)
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
    if (varLoc == "intron_variant" or varLoc == "UTR_variant") and getVarType(variant) == "substitution" and varInExon(variant) == False:
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
    Also returns sequence of closest reference splice site and genomic position of splice site
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
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"], plusStrandSeq=True)
            # splice site is 3 bp to the right of donor Start (+3 because plus strand numbering increases from left to right)
            # splice site is 3 bp to the right because exon end is 3 bp to the right of donor start
            genomicSplicePos = closestSpliceBounds["donorStart"] + 3
        else:
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"], plusStrandSeq=False)
            # splice site is 3 bp to the right of donor Start (-3 because minus strand numbering decreases from left to right)
            # splice site is 3 bp to the right because exon end is 3 bp to the right of donor start
            genomicSplicePos = closestSpliceBounds["donorStart"] - 3
        exonStart = 0
        intronStart = STD_EXONIC_PORTION
    else:
        # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
        # for plus strand it is acceptorEnd - deNovoOffset because
        # the genomic position increases from left to right on the plus strand and subtraction reduces the refSeq to correct length
        # for minus strand it is acceptorEnd + deNovoOffset because
        # the genomic position decreases from left to right on the minus strand and addition reduces the refSeq to correct length
        if getVarStrand(variant) == "+":
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"], (closestSpliceBounds["acceptorEnd"] - deNovoOffset), plusStrandSeq=True)
            # splice site is 3 bp to the left of reference acceptor End (-3 because plus strand numbering increases from left to right)
            # minus deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            genomicSplicePos = closestSpliceBounds["acceptorEnd"] - 3 - deNovoOffset
        else:
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"], (closestSpliceBounds["acceptorEnd"] + deNovoOffset), plusStrandSeq=False)
            # splice site is 3 bp to the left of reference acceptor End (+3 because minus strand numbering decreases from left to right)
            # plus deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            genomicSplicePos = closestSpliceBounds["acceptorEnd"] + 3 + deNovoOffset
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
                    "zScore": "N/A",
                    "genomicSplicePos": "N/A"}
        closestMaxEntScanScore = runMaxEntScan(refSeq, donor=donor)
        closestZScore = getZScore(closestMaxEntScanScore, donor=donor)
        return {"exonName": exonName,
                "sequence": refSeq.upper(),
                "exonStart": exonStart,
                "intronStart": intronStart,
                "maxEntScanScore": closestMaxEntScanScore,
                "zScore": closestZScore,
                "genomicSplicePos": genomicSplicePos}
    else:
        return {"exonName": exonName,
                "sequence": refSeq.upper(),
                "genomicSplicePos": genomicSplicePos}

def isCIDomainInRegion(regionStart, regionEnd, boundaries, gene):
    '''
    Given a region of interest, boundaries (either enigma or priors) and gene of interest
    Determines if there is an overlap between the region of interest and a CI domain
    For minus strand gene (BRCA1) regionStart > regionEnd
    For plus strand gene (BRCA2) regionStart < regionEnd
    Returns True if there is an overlap, False otherwise
    '''
    if gene == "BRCA1":
        if regionStart < regionEnd:
            start = regionEnd
            end = regionStart
        else:
            start = regionStart
            end = regionEnd
        for domain in brca1CIDomains[boundaries].keys():
            domainStart = brca1CIDomains[boundaries][domain]["domStart"]
            domainEnd = brca1CIDomains[boundaries][domain]["domEnd"]
            overlap = range(max(end, domainEnd), min(start, domainStart) + 1)
            if len(overlap) > 0:
                return True
    elif gene == "BRCA2":
        if regionStart < regionEnd:
            start = regionStart
            end = regionEnd
        else:
            start = regionEnd
            end = regionStart
        for domain in brca2CIDomains[boundaries].keys():
            domainStart = brca2CIDomains[boundaries][domain]["domStart"]
            domainEnd = brca2CIDomains[boundaries][domain]["domEnd"]
            overlap = range(max(start, domainStart), min(end, domainEnd) + 1)
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
        # +1 is not included in the below equation for exonLength
        # because due to RefSeq numbering varExonEnd is 1 bp too long
        # varExonEnd is first intronic base (+1 position)
        # for minus strand genes
        exonLength = varExonStart - varExonEnd
    else:
        varExonStart = int(exonBounds[varExonNum]["exonStart"])
        varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
        # +1 is not included in the below equatio for exonLength
        # because due to RefSeq numbering varExonStart is 1 bp too long
        # varExonStart is last intronic base (-1 position)
        # for plus strand genes
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
            # +1 to account for all positions including newSplicePos
            # adding one to exonStart increases length by 1 bp because numbering decreases from left to right
            exonLength = varExonStart - newSplicePos + 1
        else:
            varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
            # -1 to account for all position including newSplicePos
            # subtracting one increases length by 1 bp because numbering decreases from left to right
            exonLength = newSplicePos - varExonEnd - 1
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

def isDeNovoWildTypeSplicePosDistanceDivisibleByThree(variant, exonicPortionSize, intronicPortionSize,
                                                      deNovoDonorInRefAcc=False, donor=True):
    '''
    Given a variant, compares de novo splicing position with wild-type splicing position
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    intronicPortionSize referes to length in bp that is considered to be in intronic portion of splice site
    deNovoDonorInRefAcc argument=True if looking for de novo donor in reference splice acceptor region, False otherwise
    Donor argument determines if function is used for de novo donor (donor=True) or de novo acceptor (donor=False)
    If distance between de novo and wild-type donors is divisible by 3, returns True
    returns False otherwise
    This function is another way to check if a de novo splice site would cause a frameshift mutation
       If it returns True, then de novo splicing would not cause a frameshift
       If it returns False, then de novo splicing would cause a frameshift
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
            # +1 for minus strand donor because splice position is the right of actual splice position
            distanceBetween = deNovoSplicePos - (wildTypeSplicePos + 1)
    else:
        wildTypeSplicePos = refExonBounds[varExonNum]["exonStart"]
        if varStrand == "+":
            distanceBetween = abs(deNovoSplicePos - wildTypeSplicePos)
        else:
            # +1 for minus strand acceptor because splice position is the left of actual splice position
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
        return {"priorProb": PATHOGENIC_PROBABILITY,
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
            priorProb = PATHOGENIC_PROBABILITY
            inExonicPortionFlag = 1
            spliceRescue = 0
        else:
            inFrame = isSplicingWindowInFrame(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                              deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=True)
            # if variant causes a frameshift, no splice rescue
            if inFrame == False:
                priorProb = PATHOGENIC_PROBABILITY
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
                isDivisible = isDeNovoWildTypeSplicePosDistanceDivisibleByThree(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                                                deNovoDonorInRefAcc=deNovoDonorInRefAcc, donor=True)
                # if truncated region includes a clinically important domain or causes a frameshift
                if CIDomainInRegion == True:
                    priorProb = PATHOGENIC_PROBABILITY
                    spliceRescue = 0
                    CIDomainInRegionFlag = 1
                    inExonicPortionFlag = 0
                    frameshiftFlag = 0
                elif isDivisible == False:
                    priorProb = PATHOGENIC_PROBABILITY
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
                        priorProb = PATHOGENIC_PROBABILITY
                        spliceRescue = 0
                        lowMESFlag = 1
                    else:
                        # variant creates a stronger splice site than reference sequence
                        if altMES < LOW_MES_CUTOFF:
                            # still a weak splice site, so no change in splicing
                            priorProb = PATHOGENIC_PROBABILITY
                            spliceRescue = 0
                            lowMESFlag = 1
                        elif (altMES >= LOW_MES_CUTOFF) and (altMES <= HIGH_MES_CUTOFF):
                            # still a weak splice site, but possibility of splice rescue
                            deNovoOffset = 0
                            closestRefData = getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False,
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
                                priorProb = PATHOGENIC_PROBABILITY
                                spliceRescue = 0
                                lowMESFlag = 1
                        else:
                            # altMES > HIGH_MES_CUTOFF, strong splice site with higher possibility of splicing rescue
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
    isDivisible = isDeNovoWildTypeSplicePosDistanceDivisibleByThree(variant, STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
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
            priorProb = LOW_SPLICING_PROBABILITY
        elif (refZScore < REF_DONOR_REFZ_CUTOFF) and ((refZScore - altZScore) > MIN_REF_ALT_ZDIFF):
            priorProb = HIGH_SPLICING_PROBABILITY
        elif (refZScore < REF_DONOR_REFZ_CUTOFF) and ((refZScore - altZScore) < MIN_REF_ALT_ZDIFF):
            priorProb = MODERATE_SPLICING_PROBABILITY
        else:
            if altZScore > REF_DONOR_HIGH_CUTOFF:
                priorProb = LOW_SPLICING_PROBABILITY
            elif altZScore <= REF_DONOR_HIGH_CUTOFF and altZScore >= REF_DONOR_LOW_CUTOFF:
                priorProb = MODERATE_SPLICING_PROBABILITY
            else:
                priorProb = HIGH_SPLICING_PROBABILITY

        # capped splicing probability due to special cases of in-frame exon skipping
        varGene = variant["Gene_Symbol"]
        exonName = spliceDonorBounds["exonName"]
        if varGene == "BRCA1" and (exonName == "exon9" or exonName == "exon10"):
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY
        if varGene == "BRCA2" and exonName == "exon12":
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY
                
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
            priorProb = LOW_SPLICING_PROBABILITY
        elif (refZScore < REF_ACC_REFZ_CUTOFF) and ((refZScore - altZScore) > MIN_REF_ALT_ZDIFF):
            priorProb = HIGH_SPLICING_PROBABILITY
        elif (refZScore < REF_ACC_REFZ_CUTOFF) and ((refZScore - altZScore) < MIN_REF_ALT_ZDIFF):
            priorProb = MODERATE_SPLICING_PROBABILITY
        else:
            if altZScore > REF_ACC_HIGH_CUTOFF:
                priorProb = LOW_SPLICING_PROBABILITY
            elif altZScore <= REF_ACC_HIGH_CUTOFF and altZScore >= REF_ACC_LOW_CUTOFF:
                priorProb = MODERATE_SPLICING_PROBABILITY
            else:
                priorProb = HIGH_SPLICING_PROBABILITY

        # capped splicing probability due to special cases of in-frame exon skipping
        varGene = variant["Gene_Symbol"]
        exonName = spliceAcceptorBounds["exonName"]
        if varGene == "BRCA1" and (exonName == "exon9" or exonName == "exon10"):
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY
        if varGene == "BRCA2" and exonName == "exon12":
            if priorProb == HIGH_SPLICING_PROBABILITY:
                priorProb = CAPPED_PROBABILITY
                
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

def getDeNovoFrameshiftAndCIStatus(variant, boundaries, donor=True, deNovoDonorInRefAcc=False):
    '''
    Given a variant, boundaries (enigma or priors), donor argument, and deNovoDonorInRefAcc argument:
      donor argument = True for de novo donors, False for de novo acceptors
      deNovoDonorInRefAcc argument = True if lookign for de novo donor in ref acceptor site, False otherwise
    Determines if new splice position causes a frameshift and would disrupt a CI Domain
    If de novo splicing would cause a frameshift, returns False
    Else, checks to see if new splice position would splice out (skip) a CI domain
    If variant de novo splice position does not cause a frameshift and does not disrupt a CI domain, reutrns True
      Returns False otherwise
    '''
    frameshiftStatus = getDeNovoSpliceFrameshiftStatus(variant, donor=donor, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    # checks to make sure that variant does not cause a frameshift
    if frameshiftStatus == True:
        return False
    else:
        # determine if CI domain is in region that would be skipped by new splicing
        if varInExon(variant) == True:
            varExonNum = getVarExonNumberSNS(variant)
        else:
            if varInSpliceRegion(variant, donor=donor, deNovo=False) == True:
                spliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=False)
                varExonNum = spliceBounds["exonName"]
            else:
                if donor == True:
                    # if a variant is in an intron de novo donor cannot splice out any of the exon
                    # so no part of a CI domain will be spliced out
                    return True
        # varExonNum is a string in the format "exonN"
        varWindowPos = getVarWindowPosition(variant, donor=donor, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
        inExonicPortion = varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=donor,
                                             deNovoDonorInRefAcc=deNovoDonorInRefAcc)
        regionStart = getNewSplicePosition(variant["Pos"], getVarStrand(variant), varWindowPos, inExonicPortion,
                                           STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH, donor=donor)
        if donor == True:
            # nextExonNum parses out N from varExonNum and adds 1 to get next exon number "exonN+1"
            # uses [4:] to remove "exon" from "exonN" so can add 1 to N to get N+1
            nextExonNum = "exon" + str(int(varExonNum[4:]) + 1)
            # skips to exon 5 for any variants in BRCA1 exon 3 because exon 4 does not exist in BRCA1 RefSeq transcript
            if variant["Gene_Symbol"] == "BRCA1" and nextExonNum == "exon4":
                nextExonNum = "exon5"
            refSpliceAccBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
            regionEnd = refSpliceAccBounds[nextExonNum]["acceptorStart"]
        else:
            # prevExonNum parses out N from varExonNum and adds 1 to get previous exon number "exonN-1"
            # uses [4:] to remove "exon" from "exonN" so can subtract 1 to N to get N-1
            prevExonNum = "exon" + str(int(varExonNum[4:]) - 1)
            if variant["Gene_Symbol"] == "BRCA1" and prevExonNum == "exon4":
                prevExonNum = "exon3"
            refSpliceDonorBounds = getRefSpliceDonorBoundaries(variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
            regionEnd = refSpliceDonorBounds[prevExonNum]["donorEnd"]
        CIDomainInRegion = isCIDomainInRegion(regionStart, regionEnd, boundaries, variant["Gene_Symbol"])
        if CIDomainInRegion == False:
            return True
        return False
    
def convertGenomicPosToTranscriptPos(genomicPos, chrom, genome, transcript):
    '''
    Given a genomic position, chrom (in format "chrN"), genome (SequenceFileDB for genome),
      and transcript (pyhgvs transcript object):
    Returns a string of the transcript position at the given genomic position
    '''
    # use "T" and "A" for ref and alt because transcript position is not dependent on these values
    # converts genomic position to transcript position
    hgvs_name = str(pyhgvs.format_hgvs_name(chrom, genomicPos, "T", "A", genome, transcript))
    # parses out transcript position from full hgvs_name
    transcriptPos = str(pyhgvs.HGVSName(hgvs_name).cdna_start)
    return transcriptPos

def formatSplicePosition(position, transcript=False):
    '''
    Given a position and transcript argument, returns a formatted splice position
    If transcript = True, returns transcript formatted position "c.N"
    If transcript = False, returns genomic formatted position "g.N"
    '''
    if transcript == True:
        return "c." + str(position)
    else:
        return "g." + str(position)

def getPriorProbDeNovoDonorSNS(variant, boundaries, exonicPortionSize, genome, transcript, deNovoDonorInRefAcc=False):
    '''
    Given a variant, boundaries (either priors or enigma), and exonicPortionSize
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in an exon or is in a reference splice donor region
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing: 
      prior probability of pathogenecity and predicted qualitative engima class 
      deNovo donor MaxEntScan scores, zscores, and sequences for ref and alt
      deNovo donor genomic and transcript splice positions in format "g.N" and "c.N" respectively
      closest donor MaxEntScan score, zscore, and sequence
      closest donor genomic and transcript splice positions in format "g.N" and "c.N" respectively
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
                        "genomicSplicePos": "N/A",
                        "transcriptSplicePos": "N/A",
                        "closestGenomicSplicePos": "N/A",
                        "closestTranscriptSplicePos": "N/A",
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
            # +-1 because per HCI PRIORS website donor position is defined as being first nucleotide that is NOT included in spliced exon
            if getVarStrand(variant) == "-":
                newGenomicSplicePos = getNewSplicePosition(variant["Pos"], "-", slidingWindowInfo["varWindowPosition"],
                                                           slidingWindowInfo["inExonicPortion"], STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                            donor=True) - 1
            else:
                newGenomicSplicePos = getNewSplicePosition(variant["Pos"], "+", slidingWindowInfo["varWindowPosition"],
                                                            slidingWindowInfo["inExonicPortion"], STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                            donor=True) + 1
            deNovoOffset = 0
            subDonorInfo = getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False,
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
                priorProb = LOW_PROBABILITY
            else:
                altGreaterRefFlag = 1
                if altZScore < DE_NOVO_DONOR_LOW_CUTOFF:
                    priorProb = LOW_PROBABILITY
                elif altZScore >= DE_NOVO_DONOR_LOW_CUTOFF and altZScore < DE_NOVO_DONOR_HIGH_CUTOFF:
                    priorProb = MODERATE_DENOVO_PROBABILITY
                else:
                    priorProb = HIGH_PROBABILITY
            if (altZScore > subZScore and refAltZScore == "N/A") or (altZScore > refAltZScore and refAltZScore != "N/A"):
                # promote prior prob by one step
                if priorProb == LOW_PROBABILITY:
                    priorProb = MODERATE_DENOVO_PROBABILITY
                elif priorProb == MODERATE_DENOVO_PROBABILITY:
                    priorProb = HIGH_PROBABILITY
                else:
                    priorProb = priorProb

            if altZScore > subZScore:
                altGreaterClosestRefFlag = 1
            if altZScore > refAltZScore and refAltZScore != "N/A":
                altGreaterClosestAltFlag = 1

            if frameshiftFlag == 0 and priorProb != 0:
                frameshiftAndCIStatus = getDeNovoFrameshiftAndCIStatus(variant, boundaries, donor=True,
                                                                       deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                if frameshiftAndCIStatus == True:
                    priorProb = LOW_PROBABILITY

            # converts genomic splice position to transcript splice position
            newTranscriptSplicePos = convertGenomicPosToTranscriptPos(newGenomicSplicePos, getVarChrom(variant), genome, transcript)
            # converts closest genomic splice position to transcript splice position
            closestGenomicSplicePos = subDonorInfo["genomicSplicePos"]
            closestTranscriptSplicePos = convertGenomicPosToTranscriptPos(closestGenomicSplicePos, getVarChrom(variant), genome, transcript)
                
            return {"priorProb": priorProb,
                    "enigmaClass": getEnigmaClass(priorProb),
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
                    "genomicSplicePos": formatSplicePosition(newGenomicSplicePos, transcript=False),
                    "transcriptSplicePos": formatSplicePosition(newTranscriptSplicePos, transcript=True),
                    "closestGenomicSplicePos": formatSplicePosition(closestGenomicSplicePos, transcript=False),
                    "closestTranscriptSplicePos": formatSplicePosition(closestTranscriptSplicePos, transcript=True),
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

def getPriorProbDeNovoAcceptorSNS(variant, exonicPortionSize, deNovoLength, genome, transcript):
    '''
    Given a variant, exonic portion size, and de novo length:
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in de novo splice acceptor region 
         de novo splice acceptor region defined by deNovoLength
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing: 
      prior probability of pathogenecity and predicted qualitative engima class (both N/A)
      deNovo acceptor MaxEntScan scores, zscores, and sequence for ref and alt
      de novo acceptor genomic and transcript splice positions in format "g.N" and "c.N" respectively
      MaxEntScan score, zscore, and sequence for closest ref splice acceptor
      closest ref acceptor genomic and transcript splice positions in format "g.N" and "c.N" respectively
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
            newGenomicSplicePos = getNewSplicePosition(variant["Pos"], getVarStrand(variant), slidingWindowInfo["varWindowPosition"],
                                                       slidingWindowInfo["inExonicPortion"], STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                       donor=False)
            deNovoOffset = deNovoLength - exonicPortionSize
            closestAccInfo = getClosestSpliceSiteScores(variant, STD_DE_NOVO_OFFSET, donor=False, deNovo=True,
                                                        deNovoDonorInRefAcc=False)
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

            # converts genomic splice position to transcript splice position
            newTranscriptSplicePos = convertGenomicPosToTranscriptPos(newGenomicSplicePos, getVarChrom(variant), genome, transcript)
            # converts closest genomic splice position to transcript splice position
            closestGenomicSplicePos = closestAccInfo["genomicSplicePos"]
            closestTranscriptSplicePos = convertGenomicPosToTranscriptPos(closestGenomicSplicePos, getVarChrom(variant), genome, transcript)
                                                                          
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
                    "genomicSplicePos": formatSplicePosition(newGenomicSplicePos, transcript=False),
                    "transcriptSplicePos": formatSplicePosition(newTranscriptSplicePos, transcript=True),
                    "closestGenomicSplicePos": formatSplicePosition(closestGenomicSplicePos, transcript=False),
                    "closestTranscriptSplicePos": formatSplicePosition(closestTranscriptSplicePos, transcript=True),
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

def getPriorProbSpliceDonorSNS(variant, boundaries, variantData, genome, transcript):
    '''
    Given a variant, boundaries (either PRIORS or ENIGMA), and a list of dictionaries with variant data
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
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
        deNovoSpliceInfo = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                      deNovoDonorInRefAcc=False)
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
                "deNovoDonorGenomicSplicePos": deNovoSpliceInfo["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoSpliceInfo["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoSpliceInfo["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoSpliceInfo["closestTranscriptSplicePos"],
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
                "deNovoAccGenomicSplicePos": "N/A",
                "deNovoAccTranscriptSplicePos": "N/A",
                "closestAccGenomicSplicePos": "N/A",
                "closestAccTranscriptSplicePos": "N/A",
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

def getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData, genome, transcript):
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
        deNovoAccInfo = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, genome, transcript)
        refPrior = refSpliceInfo["priorProb"]
        proteinPrior = "N/A"
        applicablePrior = refSpliceInfo["priorProb"]
        if varInExon(variant) == True:
            deNovoDonorInfo = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                         deNovoDonorInRefAcc=True)
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
                               "genomicSplicePos": "N/A",
                               "transcriptSplicePos": "N/A",
                               "closestGenomicSplicePos": "N/A",
                               "closestTranscriptSplicePos": "N/A",
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
                "deNovoDonorGenomicSplicePos": deNovoDonorInfo["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorInfo["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorInfo["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorInfo["closestTranscriptSplicePos"],
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
                "deNovoAccGenomicSplicePos": deNovoAccInfo["genomicSplicePos"],
                "deNovoAccTranscriptSplicePos": deNovoAccInfo["transcriptSplicePos"],
                "closestAccGenomicSplicePos": deNovoAccInfo["closestGenomicSplicePos"],
                "closestAccTranscriptSplicePos": deNovoAccInfo["closestTranscriptSplicePos"],
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
        varHGVS = variant["HGVS_cDNA"]
        if varHGVS == "-":
            # if HGVS_cDNA field is blank use pyhgvs field instead
            # [12:] parses out the NM_ accession so varHGVS is in format c.65C>T
            varHGVS = variant["pyhgvs_cDNA"][12:]
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
        if proteinPrior == PATHOGENIC_PROBABILITY:
            proteinPrior = CAPPED_PROBABILITY

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
                "deNovoDonorGenomicSplicePos": "-",
                "deNovoDonorTranscriptSplicePos": "-",
                "closestDonorGenomicSplicePos": "-",
                "closestDonorTranscriptSplicePos": "-",
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
                "deNovoAccGenomicSplicePos": "-",
                "deNovoAccTranscriptSplicePos": "-",
                "closestAccGenomicSplicePos": "-",
                "closestAccTranscriptSplicePos": "-",
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
    
def getPriorProbInExonSNS(variant, boundaries, variantData, genome, transcript):
    '''
    Given a variant, boundaries (either "enigma" or "priors") and a list of dictionaries containing variant data:
      1. Checks that variant is in an exon or clinically important domains and NOT in a splice site
      2. Checks that variant is a SNS variant
      3. Gets protein prior from variantData
      4. Determines if variant is a nonsense variant, if yes determines if splice rescue occurs
      5. If not a nonsense variant, calculates de novo donor prior and de novo acceptor prior if applicable
         Gets applicable prior if variant has a de novo donor prior
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
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
        deNovoDonorData = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                     deNovoDonorInRefAcc=False)
        if varInSpliceRegion(variant, donor=False, deNovo=True):
            deNovoAccData = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, genome, transcript)
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
                             "genomicSplicePos": "N/A",
                             "transcriptSplicePos": "N/A",
                             "closestGenomicSplicePos": "N/A",
                             "closestTranscriptSplicePos": "N/A",
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
                "deNovoDonorGenomicSplicePos": deNovoDonorData["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorData["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorData["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorData["closestTranscriptSplicePos"],
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
                "deNovoAccGenomicSplicePos": deNovoAccData["genomicSplicePos"],
                "deNovoAccTranscriptSplicePos": deNovoAccData["transcriptSplicePos"],
                "closestAccGenomicSplicePos": deNovoAccData["closestGenomicSplicePos"],
                "closestAccTranscriptSplicePos": deNovoAccData["closestTranscriptSplicePos"],
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
        priorProb = LOW_PROBABILITY
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
                "deNovoDonorGenomicSplicePos": "-",
                "deNovoDonorTranscriptSplicePos": "-",
                "closestDonorGenomicSplicePos": "-",
                "closestDonorTranscriptSplicePos": "-",
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
                "deNovoAccGenomicSplicePos": "-",
                "deNovoAccTranscriptSplicePos": "-",
                "closestAccGenomicSplicePos": "-",
                "closestAccTranscriptSplicePos": "-",
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

def getPriorProbIntronicDeNovoDonorSNS(variant, genome, transcript):
    '''
    Given a variant,
      1. Checks that variant is NOT in exon or reference donor/acceptor site
      2. Checks that variant is a substitution variant
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Determines if alt MES score is greater than ref MES score for highest scoring sliding window
      If altMES > refMES, then deNovoDonorAltGreaterRefFlag = 1 (0 otherwise)
    Determines closest ref donor scores
      If altMES > closestRefDonorMES, then deNovoDonorAltGreaterClosestRefFlag = 1 (0 otherwise)
    If either flag is equal to 1, flag variant for further analysis (spliceFlag = 1), spliceFlag = 0 otherwise
    altGreaterClosestAltFlag is always equals N/A because this function is not used for any variants in ref splice sites
    Returns dictionary containing prior prob, enigma class, de novo donor scores, and splice flag
    Also contains closest ref donor scores and de novo and closest donor genomic and transcript splice positions
    '''
    inExon = varInExon(variant)
    inRefDonor = varInSpliceRegion(variant, donor=True, deNovo=False)
    inRefAcc = varInSpliceRegion(variant, donor=False, deNovo=False)
    if inExon == False and inRefDonor == False and inRefAcc == False:
        if getVarType(variant) == "substitution":
            deNovoDonorInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                                    donor=True, deNovo=False, deNovoDonorInRefAcc=False)
            # +-1 because per HCI PRIORS website donor position is defined as being first nucleotide that is NOT included in spliced exon
            if getVarStrand(variant) == "-":
                 newGenomicSplicePos = getNewSplicePosition(variant["Pos"], "-", deNovoDonorInfo["varWindowPosition"],
                                                            deNovoDonorInfo["inExonicPortion"], STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                            donor=True) - 1
            else:
                newGenomicSplicePos = getNewSplicePosition(variant["Pos"], "+", deNovoDonorInfo["varWindowPosition"],
                                                            deNovoDonorInfo["inExonicPortion"], STD_EXONIC_PORTION, STD_ACC_INTRONIC_LENGTH,
                                                            donor=True) + 1
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

            # converts genomic splice position to transcript splice position
            newTranscriptSplicePos = convertGenomicPosToTranscriptPos(newGenomicSplicePos, getVarChrom(variant), genome, transcript)
            # converts closest genomic splice position to transcript splice position
            closestGenomicSplicePos = closestDonorInfo["genomicSplicePos"]
            closestTranscriptSplicePos = convertGenomicPosToTranscriptPos(closestGenomicSplicePos, getVarChrom(variant), genome, transcript)

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
                    "genomicSplicePos": formatSplicePosition(newGenomicSplicePos, transcript=False),
                    "transcriptSplicePos": formatSplicePosition(newTranscriptSplicePos, transcript=True),
                    "closestGenomicSplicePos": formatSplicePosition(closestGenomicSplicePos, transcript=False),
                    "closestTranscriptSplicePos": formatSplicePosition(closestTranscriptSplicePos, transcript=True),
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
    
def getPriorProbInIntronSNS(variant, boundaries, genome, transcript):
    '''
    Given a variant and boundaries (either "priors or "enigma"),
    Checks that variant is located in an intron and is a substitution variant
    Determines if variant creates a de novo donor site in the intron
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
    Returns a dictionary containing applicable prior and predicted qualitative enigma class
    Dictionary also contains de novo donor ref and alt scores
    AND a spliceFlag which is equal to 1 if variant creates a better de novo splice site than ref sequence
    Rest of values in dictionary are equal to 0, "-", or N/A because they are not relevant to variant
    '''
    varLoc = getVarLocation(variant, boundaries)
    varType = getVarType(variant)
    if varLoc == "intron_variant" and varType == "substitution":
        deNovoDonorData = getPriorProbIntronicDeNovoDonorSNS(variant, genome, transcript)
        if deNovoDonorData["spliceFlag"] == 0:
            priorProb = LOW_PROBABILITY
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
                "deNovoDonorGenomicSplicePos": deNovoDonorData["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorData["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorData["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorData["closestTranscriptSplicePos"],
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
                "deNovoAccGenomicSplicePos": "N/A",
                "deNovoAccTranscriptSplicePos": "N/A",
                "closestAccGenomicSplicePos": "N/A",
                "closestAccTranscriptSplicePos": "N/A",
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

def getPriorProbInUTRSNS(variant, boundaries, genome, transcript):
    '''
    Given a variant and boundaries (either "priors" or "enigma"),
    Checks that variant is a SNS variant in a UTR
    Determines prior prob based on location (5'/3' UTR and intron/exon)
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
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
                         "genomicSplicePos": "N/A",
                         "transcriptSplicePos": "N/A",
                         "closestGenomicSplicePos": "N/A",
                         "closestTranscriptSplicePos": "N/A",
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
                           "genomicSplicePos": "N/A",
                           "transcriptSplicePos": "N/A",
                           "closestGenomicSplicePos": "N/A",
                           "closestTranscriptSplicePos": "N/A",
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
            applicablePrior = LOW_PROBABILITY
            applicableClass = getEnigmaClass(applicablePrior)
            spliceFlag = 0
        elif varCons == "5_prime_UTR_variant":
            if varInExon(variant) == True:
                deNovoDonorData = getPriorProbDeNovoDonorSNS(variant, boundaries, STD_EXONIC_PORTION, genome, transcript,
                                                             deNovoDonorInRefAcc=False)
                applicablePrior = deNovoDonorData["priorProb"]
                applicableClass = deNovoDonorData["enigmaClass"]
                spliceFlag = 0
                if varInSpliceRegion(variant, donor=False, deNovo=True) == True:
                    deNovoAccData = getPriorProbDeNovoAcceptorSNS(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                                  genome, transcript)
        else:
            # to account for variants in 5' UTR that are classified as other variant types
            deNovoDonorData = getPriorProbIntronicDeNovoDonorSNS(variant, genome, transcript)
            spliceFlag = deNovoDonorData["spliceFlag"]
            if spliceFlag == 1:
                applicablePrior = deNovoDonorData["priorProb"]
                applicableClass = deNovoDonorData["enigmaClass"]
            else:
                applicablePrior = LOW_PROBABILITY
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
                "deNovoDonorGenomicSplicePos": deNovoDonorData["genomicSplicePos"],
                "deNovoDonorTranscriptSplicePos": deNovoDonorData["transcriptSplicePos"],
                "closestDonorGenomicSplicePos": deNovoDonorData["closestGenomicSplicePos"],
                "closestDonorTranscriptSplicePos": deNovoDonorData["closestTranscriptSplicePos"],
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
                "deNovoAccGenomicSplicePos": deNovoAccData["genomicSplicePos"],
                "deNovoAccTranscriptSplicePos": deNovoAccData["transcriptSplicePos"],
                "closestAccGenomicSplicePos": deNovoAccData["closestGenomicSplicePos"],
                "closestAccTranscriptSplicePos": deNovoAccData["closestTranscriptSplicePos"],
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

def getVarData(variant, boundaries, variantData, genome, transcript):
    '''
    Given variant, boundaries (either "priors" or "enigma') and list of dictionaries with variant data
    Genome is a SequenceFileDB for genome and transcript is a pyhgvs transcript object)
       both genome and transcript are necessary to convert from genomic to transcript coordinates
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
                 "deNovoDonorGenomicSplicePos": "-",
                 "deNovoDonorTranscriptSplicePos": "-",
                 "closestDonorGenomicSplicePos": "-",
                 "closestDonorTranscriptSplicePos": "-",
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
                 "deNovoAccGenomicSplicePos": "-",
                 "deNovoAccTranscriptSplicePos": "-",
                 "closestAccGenomicSplicePos": "-",
                 "closestAccTranscriptSplicePos": "-",
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
                varData =  getPriorProbSpliceDonorSNS(variant, boundaries, variantData, genome, transcript)
            elif varLoc == "CI_splice_acceptor_variant" or varLoc == "splice_acceptor_variant":
                varData =  getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData, genome, transcript)
            elif varLoc == "CI_domain_variant" or varLoc == "exon_variant":
                varData = getPriorProbInExonSNS(variant, boundaries, variantData, genome, transcript)
            elif varLoc == "grey_zone_variant":
                varData = getPriorProbInGreyZoneSNS(variant, boundaries, variantData)
            elif varLoc == "after_grey_zone_variant":
                varData =  getPriorProbAfterGreyZoneSNS(variant, boundaries)
            elif varLoc == "UTR_variant":
                varData =  getPriorProbInUTRSNS(variant, boundaries, genome, transcript)
            elif varLoc == "intron_variant":
                varData =  getPriorProbInIntronSNS(variant, boundaries, genome, transcript)
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
    parser.add_argument('-g', "--genomeFile", help="Fasta file containing hg38 reference genome")
    parser.add_argument('-t', "--transcriptFile", help="RefSeq annotation hg38-based genepred file")
    args = parser.parse_args()    

    inputData = csv.DictReader(open(args.inputFile, "r"), delimiter="\t")
    fieldnames = inputData.fieldnames
    newHeaders = ["varType", "varLoc", "applicablePrior", "applicableEnigmaClass", "proteinPrior", "refDonorPrior", "deNovoDonorPrior",
                  "refRefDonorMES", "refRefDonorZ", "altRefDonorMES", "altRefDonorZ", "refRefDonorSeq", "altRefDonorSeq", "refDonorVarStart",
                  "refDonorVarLength", "refDonorExonStart", "refDonorIntronStart", "refDeNovoDonorMES", "refDeNovoDonorZ", "altDeNovoDonorMES",
                  "altDeNovoDonorZ", "refDeNovoDonorSeq", "altDeNovoDonorSeq", "deNovoDonorVarStart", "deNovoDonorVarLength",
                  "deNovoDonorExonStart", "deNovoDonorIntronStart", "deNovoDonorGenomicSplicePos", "deNovoDonorTranscriptSplicePos",
                  "closestDonorGenomicSplicePos", "closestDonorTranscriptSplicePos", "closestDonorRefMES", "closestDonorRefZ",
                  "closestDonorRefSeq", "closestDonorAltMES", "closestDonorAltZ", "closestDonorAltSeq", "closestDonorExonStart",
                  "closestDonorIntronStart", "deNovoDonorAltGreaterRefFlag", "deNovoDonorAltGreaterClosestRefFlag",
                  "deNovoDonorAltGreaterClosestAltFlag", "deNovoDonorFrameshiftFlag", "refAccPrior", "deNovoAccPrior", "refRefAccMES",
                  "refRefAccZ", "altRefAccMES", "altRefAccZ", "refRefAccSeq", "altRefAccSeq", "refAccVarStart", "refAccVarLength", "refAccExonStart",
                  "refAccIntronStart", "refDeNovoAccMES", "refDeNovoAccZ", "altDeNovoAccMES", "altDeNovoAccZ", "refDeNovoAccSeq", "altDeNovoAccSeq",
                  "deNovoAccVarStart", "deNovoAccVarLength", "deNovoAccExonStart", "deNovoAccIntronStart", "deNovoAccGenomicSplicePos",
                  "deNovoAccTranscriptSplicePos", "closestAccGenomicSplicePos", "closestAccTranscriptSplicePos", "closestAccRefMES",
                  "closestAccRefZ", "closestAccRefSeq", "closestAccAltMES", "closestAccAltZ", "closestAccAltSeq",
                  "closestAccExonStart", "closestAccIntronStart", "deNovoAccAltGreaterRefFlag", "deNovoAccAltGreaterClosestRefFlag",
                  "deNovoAccAltGreaterClosestAltFlag", "deNovoAccFrameshiftFlag", "spliceSite", "spliceRescue", "spliceFlag", "frameshiftFlag",
                  "inExonicPortionFlag", "CIDomainInRegionFlag", "isDivisibleFlag", "lowMESFlag"]
    for header in newHeaders:
        fieldnames.append(header)
    outputData = csv.DictWriter(open(args.outputFile, "w"), delimiter="\t", fieldnames=fieldnames)
    outputData.writerow(dict((fn,fn) for fn in inputData.fieldnames))

    # read genome sequence
    genome38 = SequenceFileDB(args.genomeFile)

    # read RefSeq transcripts
    with open(args.transcriptFile) as infile:
        transcripts = pyhgvs_utils.read_transcripts(infile)

    def get_transcript(name):
        return transcripts.get(name)

    brca1Transcript = get_transcript(BRCA1_RefSeq)
    brca2Transcript = get_transcript(BRCA2_RefSeq)

    totalVariants = 0
    for variant in inputData:
        variantData = csv.DictReader(open(args.variantFile, "r"), delimiter="\t")
        if variant["Gene_Symbol"] == "BRCA1":
            varData = getVarData(variant, args.boundaries, variantData, genome38, brca1Transcript)
        elif variant["Gene_Symbol"] == "BRCA2":
            varData = getVarData(variant, args.boundaries, variantData, genome38, brca2Transcript)
        variant = addVarDataToRow(varData, variant)
        outputData.writerow(variant)
        totalVariants += 1
        print "variant", totalVariants, "complete"
    
if __name__ == "__main__":
    main()
