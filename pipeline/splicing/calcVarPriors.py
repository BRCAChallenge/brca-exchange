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

def getFastaSeq(chrom, rangeStart, rangeStop):
    '''
    Given chromosome (in format 'chr13'), region genomic start position, and
    region genomic end position:
    Returns a string containing the sequence inclusive of rangeStart and rangeStop
    To get + strand sequence, rangeStart < rangeStop
    '''
    url = "http://togows.org/api/ucsc/hg38/%s:%d-%d.fasta" % (chrom, rangeStart, rangeStop)
    req = requests.get(url)
    lines = req.content.split('\n')

    sequence = ""
    for base in range(1, len(lines)):
        sequence += lines[base]
    return sequence

def getSeqLocDict(chrom, strand, rangeStart, rangeStop):
    '''
    Given chromosome, strand, region genomic start position, and region genomic end position
    returns a dictionary containing the genomic position as the key and reference allele as the value
    For BRCA1 rangeStart > rangeStop, for BRCA2 rangeStart < rangeEnd
    '''
    seqLocDict = {}
    if strand == "-":
        regionStart = int(rangeStop)
        regionEnd = int(rangeStart)
    else:
        regionStart = int(rangeStart)
        regionEnd = int(rangeStop)

    sequence = getFastaSeq(chrom, regionStart, regionEnd)
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
        regionStart = rangeStop
        regionEnd = rangeStart
    else:
        regionStart = rangeStart
        regionEnd = rangeStop
    refSeq = getFastaSeq(varChrom, regionStart, regionEnd)
    if varStrand == "-":
        refSeq = str(Seq(refSeq).reverse_complement())
    refSeqDict = getSeqLocDict(varChrom, varStrand, rangeStart, rangeStop)
    altSeqDict = getAltSeqDict(variant, refSeqDict)
    altSeq = getAltSeq(altSeqDict, varStrand)
    return {"refSeq": refSeq,
            "altSeq": altSeq}    
    
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
            refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["donorStart"], refSpliceBounds["donorEnd"]).upper()
        else:
            refSpliceBounds = getVarSpliceRegionBounds(variant, donor=donor, deNovo=True)
            deNovoOffset = deNovoLength - exonicPortionSize
            # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
            if getVarStrand(variant) == "+":
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] - deNovoOffset)).upper()
            else:
                refSpliceSeq = getFastaSeq(getVarChrom(variant), refSpliceBounds["acceptorStart"],
                                           (refSpliceBounds["acceptorEnd"] + deNovoOffset)).upper()
        for position, seqs in slidingWindowInfo["windowSeqs"].iteritems():
            if seqs["refSeq"] == refSpliceSeq:
                refSpliceWindow = position
                # removes reference splice window so it is not considered for de novo splicing
                del windowAltMaxEntScanScores[refSpliceWindow]
    # to get tuple containing sequence with variant position with maximum alt MaxEntScan score
    maxAltWindowScore = max(windowAltMaxEntScanScores.items(), key=lambda k: k[1])
    maxVarPosition = maxAltWindowScore[0]
    maxScores = slidingWindowInfo["windowScores"][maxVarPosition]
    
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

def getClosestSpliceSiteScores(variant, deNovoOffset, donor=True, deNovo=False, deNovoDonorInRefAcc=False, testMode=False):
    '''
    Given a variant, determines scores for closest reference splice sequence
    deNovoOffset refers to difference between de novo acceptor length and exonic portion size
       If donor = True, looks for closest splice donor sequence
       If donor = False, looks for closest splice acceptor sequence
       If deNovo = True, accomodates for de novo splicing
         *Note only use argument deNovo=True in this function if donor=False
         *Function will not return correct sequence if donor=True and deNovo=True
    If exonic variant, returns a dictionary containing:
       MaxEntScan score and z-score for reference closest splice sequence
    If variant located in referene splice site, returns a dictionary containing:
       MaxEntScan score and z-score for that reference splice site sequence
    deNovoDonorInRefAcc = False if NOT checking for de novo splice donor sites in reference splice acceptor sites
    deNovoDonorInRefAcc = True if checking for de novo splice donor sites in reference splice acceptor sites 
    '''
    varGenPos = int(variant["Pos"])
    varChrom = getVarChrom(variant)
    if varInExon(variant) == True and deNovo == False:
        exonNumber = getVarExonNumberSNS(variant)
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
        refSeq = getFastaSeq(varChrom, closestSpliceBounds["donorStart"], closestSpliceBounds["donorEnd"])
    else:
        # acceptorEnd +- deNovoOffset because deNovo splice acceptor region is deNovoOffset bp longer than reference splice acceptor region
        if getVarStrand(variant) == "+":
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"], (closestSpliceBounds["acceptorEnd"] - deNovoOffset))
        else:
            refSeq = getFastaSeq(varChrom, closestSpliceBounds["acceptorStart"], (closestSpliceBounds["acceptorEnd"] + deNovoOffset))
    if testMode == False:
        # to prevent issue with running max ent scan score on unittests
        closestMaxEntScanScore = runMaxEntScan(refSeq, donor=donor)
        closestZScore = getZScore(closestMaxEntScanScore, donor=donor)
        return {"exonName": exonName,
                "maxEntScanScore": closestMaxEntScanScore,
                "zScore": closestZScore}
    else:
        return {"exonName": exonName}

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

def getRefExonLength(variant):
    '''Given a variant, returns the length of the reference exon'''
    if varInExon(variant) == True:
        varExonNum = getVarExonNumberSNS(variant)
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

def getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, exonicPortionSize):
    '''
    Given a variant's:
    genetic postion, strand, sliding window position with max MES score AND
      whether that position is within exonic portion of highest scoring window and exonic portion size
    Returns the position where splicing occurs for a de novo splice donor
    '''
    if varStrand == "+":
        if inExonicPortion == False:
            newSplicePos = int(varGenPos) - (varWindowPos - exonicPortionSize)
        else:
            newSplicePos = int(varGenPos) + abs(varWindowPos - exonicPortionSize)
    else:
        if inExonicPortion == False:
            newSplicePos = int(varGenPos) + (varWindowPos - exonicPortionSize)
        else:
            newSplicePos = int(varGenPos) - abs(varWindowPos - exonicPortionSize)
    return newSplicePos
    
def getAltExonLength(variant, exonicPortionSize, deNovoDonorInRefAcc=False):
    '''
    Given a variant and the exonic portion size,
    returns the length of the alternate exon after splicing occurs in max MES window
    Function can only be used for de novo donor variants
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    '''
    if varInExon(variant) == True:
        varExonNum = getVarExonNumberSNS(variant)
        exonBounds = getExonBoundaries(variant)
        slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH, donor=True,
                                                                  deNovoDonorInRefAcc=deNovoDonorInRefAcc)
        newSplicePos = getNewSplicePosition(variant["Pos"], getVarStrand(variant), slidingWindowInfo["varWindowPosition"],
                                            slidingWindowInfo["inExonicPortion"], exonicPortionSize)
        if getVarStrand(variant) == "-":
            varExonStart = int(exonBounds[varExonNum]["exonStart"])
            # newSplicePos -1 to account for RefSeq numbering which starts to the right of the first base
            exonLength = varExonStart - (newSplicePos - 1)
        else:
            varExonEnd = int(exonBounds[varExonNum]["exonEnd"])
            refExonLength = getRefExonLength(variant)
            # need to compare to refExonLength because of + strand gene
            exonLength = refExonLength - (varExonEnd - newSplicePos)
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

def isSplicingWindowInFrame(variant, exonicPortionSize, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determines ref and alt exon length and compares them
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    deNovoDonorInRefAcc=True if looking for deNovoDonor in ref acceptor site, False otherwise
    If ref and alt exon are in the same reading frame, returns True
    '''
    refLength = getRefExonLength(variant)
    altLength = getAltExonLength(variant, exonicPortionSize, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
    return compareRefAltExonLengths(refLength, altLength)

def compareDeNovoWildTypeSplicePos(variant, exonicPortionSize, deNovoDonorInRefAcc=False):
    '''
    Given a variant, compares de novo splicing position with wild-type splicing position
    exonicPortionSize refers to length in bp that is considered to be in exonic portion of splice site
    deNovoDonorInRefAcc argument=True if looking for de novo donor in reference splice acceptor region, False otherwise
    If distance between de novo and wild-type donors is divisible by 3, returns True
    returns False otherwise
    '''
    if varInExon(variant) == True:
        varStrand = getVarStrand(variant)
        varExonNum = getVarExonNumberSNS(variant)
        refExonBounds = getExonBoundaries(variant)
        wildTypeSplicePos = refExonBounds[varExonNum]["exonEnd"]
        slidingWindowInfo = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH, donor=True,
                                                                  deNovoDonorInRefAcc=deNovoDonorInRefAcc)
        deNovoSplicePos = getNewSplicePosition(variant["Pos"], varStrand, slidingWindowInfo["varWindowPosition"],
                                               slidingWindowInfo["inExonicPortion"], exonicPortionSize)
        if varStrand == "+":
            distanceBetween = wildTypeSplicePos - deNovoSplicePos
        else:
            # +1 because of RefSeq numbering
            distanceBetween = deNovoSplicePos - (wildTypeSplicePos + 1)
        if distanceBetween % 3 == 0:
            return True
        else:
            return False
        
def getPriorProbSpliceRescueNonsenseSNS(variant, boundaries, deNovoDonorInRefAcc=False):
    '''
    Given a variant, determines if there is a possibility of splice rescue
    deNovoDonorInRefAcc argument = True  if looking for deNovoDonor in ref acceptor site, False otherwise
    If there is a possibility of splice rescue, flags variant for further analysis
    Else assigns prior probability of pathogenecity and predicted qualitative ENIGMA class
    '''
    varCons = getVarConsequences(variant)
    # check that variant causes a premature stop codon in an exon
    if varCons == "stop_gained" and varInExon(variant) == True:
        spliceFlag = 0
        spliceRescue = 0
        frameshift = 0
        # if variant is in specified exonic portion of highest scoring sliding window, no splice rescue
        if varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=True,
                              deNovoDonorInRefAcc=deNovoDonorInRefAcc) == True:
            priorProb = 0.97
            spliceRescue = 0
            frameshift = 0
        else:
            inFrame = isSplicingWindowInFrame(variant, STD_EXONIC_PORTION, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            # if variant causes a frameshift, no splice rescue
            if inFrame == False:
                priorProb = 0.99
                spliceRescue = 0
                frameshift = 1
            else:
                varGenPos = int(variant["Pos"])
                varStrand = getVarStrand(variant)
                varExonNum = getVarExonNumberSNS(variant)
                # varExonNum returns a string in the format "exonN"
                # nextExonNum parses out N from varExonNum and adds 1 to get next exon number key "exonN+1"
                # use [4:] to remove "exon" from "exonN" so can add 1 to N to get N+1
                nextExonNum = "exon" + str(int(varExonNum[4:]) + 1)
                refSpliceAccBounds = getSpliceAcceptorBoundaries(variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
                varWindowPos = getVarWindowPosition(variant, donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                inExonicPortion = varInExonicPortion(variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, donor=True,
                                                     deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                # gets region from new splice position to next splice acceptor
                regionStart = getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, STD_EXONIC_PORTION)
                regionEnd = refSpliceAccBounds[nextExonNum]["acceptorStart"]
                CIDomainInRegion = isCIDomainInRegion(regionStart, regionEnd, boundaries, variant["Gene_Symbol"])
                isDivisible = compareDeNovoWildTypeSplicePos(variant, STD_EXONIC_PORTION, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
                # if truncated region includes a clinically important domain or causes a frameshift
                if CIDomainInRegion == True or isDivisible == False:
                    priorProb = 0.99
                    spliceRescue = 0
                    if isDivisible == False:
                        frameshift = 1
                else:
                    # possibility of splice rescue, flag variant for further splicing assays
                    spliceFlag = 1
                    spliceRescue = 1
                    frameshift = 0
                    priorProb = "N/A"
                    enigmaClass = "N/A"

        if priorProb == "N/A":
            pass
        else:
            enigmaClass = getEnigmaClass(priorProb)

        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass,
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshift": frameshift}

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
     prior probability of pathogenecity, predicted qualitative enigma class, and ref and alt MES and zscores
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
     prior probability of pathogenecity, predicted qualitative enigma class, and ref and alt MES and zscores
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
    '''
    varType = getVarType(variant)
    varLoc = getVarLocation(variant, boundaries)
    varCons = getVarConsequences(variant)
    if varType == "substitution" and varLoc == "after_grey_zone_variant":
        if varCons == "stop_gained" or varCons == "missense_variant":
            priorProb = "N/A"
            enigmaClass = "class_2"
        return {"priorProb": priorProb,
                "enigmaClass": enigmaClass}

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
      deNovo donor MaxEntScan scores and zscores for ref and alt
    deNovoDonorInRefAcc = False if NOT checking for de novo donor in reference splice acceptor site
    deNovoDonorInRefAcc = True if checking for de novo donors in reference splice acceptor site
    '''
    if getVarType(variant) == "substitution":
        if varInSpliceRegion(variant, donor=True, deNovo=True) == True:
            if varInExon(variant) == True and varInIneligibleDeNovoExon(variant, donor=True) == True:
                return {"priorProb": "N/A",
                        "enigmaClass": "N/A",
                        "refMaxEntScanScore": "-",
                        "altMaxEntScanScore": "-",
                        "refZScore": "-",
                        "altZScore": "-",
                        "deNovoDonorFlag": 0}
            slidingWindowScores = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, STD_DE_NOVO_LENGTH,
                                                                        donor=True, deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            subDonorScores = getClosestSpliceSiteScores(variant, STD_DE_NOVO_OFFSET, donor=True, deNovo=False,
                                                        deNovoDonorInRefAcc=deNovoDonorInRefAcc)
            altZScore = slidingWindowScores["altZScore"]
            refZScore = slidingWindowScores["refZScore"]
            if altZScore <= refZScore:
                priorProb = 0
            else:
                if altZScore < -2.0:
                    priorProb = 0.02
                elif altZScore >= -2.0 and altZScore < 0.0:
                    priorProb = 0.3
                else:
                    priorProb = 0.64
            if altZScore > subDonorScores["zScore"]:
                # promote prior prob by one step
                if priorProb == 0:
                    priorProb = 0.3
                elif priorProb == 0.02:
                    priorProb = 0.3
                elif priorProb == 0.3:
                    priorProb = 0.64
                else:
                    priorProb = priorProb
            
            if priorProb == 0: 
                priorProb = "N/A"
                enigmaClass = "N/A"
                deNovoDonorFlag = 0
            else:
                priorProb = priorProb
                enigmaClass = getEnigmaClass(priorProb)
                deNovoDonorFlag = 1
                
            return {"priorProb": priorProb,
                    "enigmaClass": enigmaClass,
                    "refMaxEntScanScore": slidingWindowScores["refMaxEntScanScore"],
                    "altMaxEntScanScore": slidingWindowScores["altMaxEntScanScore"],
                    "refZScore": refZScore,
                    "altZScore": altZScore,
                    "deNovoDonorFlag": deNovoDonorFlag}

def getPriorProbDeNovoAcceptorSNS(variant, exonicPortionSize, deNovoLength):
    '''
    Given a variant, exonic portion size, and de novo length:
      1. checks that variant is a single nucleotide substitution
      2. checks that variant is in de novo splice acceptor region 
         de novo splice acceptor region defined by deNovoLength
    Returns a dictionary containing: 
      prior probability of pathogenecity and predicted qualitative engima class (both N/A)
      deNovo acceptor MaxEntScan scores and zscores for ref and alt sequence
      deNovo acceptor flag: which equals 1 alt > ref
    '''
    if getVarType(variant) == "substitution":
        if varInSpliceRegion(variant, donor=False, deNovo=True) == True:
            slidingWindowScores = getMaxMaxEntScanScoreSlidingWindowSNS(variant, exonicPortionSize, deNovoLength,
                                                                        donor=False)
            deNovoOffset = deNovoLength - exonicPortionSize
            closestAccScores = getClosestSpliceSiteScores(variant, deNovoOffset, donor=False, deNovo=True)
            altZScore = slidingWindowScores["altZScore"]
            refZScore = slidingWindowScores["refZScore"]
            closestZScore = closestAccScores["zScore"]
            if altZScore <= refZScore:
                deNovoAccFlag = 0
            else:
                deNovoAccFlag = 1
                        
            return {"priorProb": "N/A",
                    "enigmaClass": "N/A",
                    "refMaxEntScanScore": slidingWindowScores["refMaxEntScanScore"],
                    "altMaxEntScanScore": slidingWindowScores["altMaxEntScanScore"],
                    "refZScore": refZScore,
                    "altZScore": altZScore,
                    "deNovoAccFlag": deNovoAccFlag}

def getPriorProbSpliceDonorSNS(variant, boundaries, variantData):
    '''
    Given a variant, boundaries (either PRIORS or ENIGMA), and a list of dictionaries with variant data
    Determines reference donor and de novo donor scores for variant
    If variant causes a nonsense mutation, determines if splice rescue occurs
    Returns dicitionary containing scores for ref and de novo splice donor/acceptor
        and protein prior if variant in exon
        score = "-" if score not applicable for variant
    Also contains other values:
        applicable prior, highest prior if variant has multiple priors
        applicable prior, highest prior if variant has multiple priors
        ref prior, prior prob for reference splice donor
        de novo prior, prior prob for de novo donor sequence
        splice flag = 1, because variant in reference splice site
        de novo donor flag = 1 if variant is possible de novo donor variant
        de novo acc flag = 0, because not applicable for variants in ref splice sites
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshift = 1 if nonsense variant causes a frameshift mutation also, 0 otherwise
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

        spliceRescue = 0
        spliceFlag = 0
        frameshift = 0
        # to check for nonsense variants in exonic portion of splice donor site
        if varInExon(variant) == True and getVarConsequences(variant) == "stop_gained":
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries)
            applicablePrior = nonsenseData["priorProb"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshift = nonsenseData["frameshift"]
            
        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": getEnigmaClass(applicablePrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": refPrior,
                "deNovoDonorPrior": deNovoPrior,
                "refRefDonorMES": refSpliceInfo["refMaxEntScanScore"],
                "refRefDonorZ": refSpliceInfo["refZScore"],
                "altRefDonorMES": refSpliceInfo["altMaxEntScanScore"],
                "altRefDonorZ": refSpliceInfo["altZScore"],
                "refDeNovoDonorMES": deNovoSpliceInfo["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoSpliceInfo["refZScore"],
                "altDeNovoDonorMES": deNovoSpliceInfo["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoSpliceInfo["altZScore"],
                "deNovoDonorFlag": deNovoSpliceInfo["deNovoDonorFlag"],
                "refAccPrior": "N/A",
                "deNovoAccPrior": "N/A",
                "refRefAccMES": "-",
                "refRefAccZ": "-",
                "altRefAccMES": "-",
                "altRefAccZ": "-",
                "refDeNovoAccMES": "-",
                "refDeNovoAccZ": "-",
                "altDeNovoAccMES": "-",
                "altDeNovoAccZ": "-",
                "deNovoAccFlag": 0,
                "spliceSite": refSpliceInfo["spliceSite"],
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshift": frameshift}

def getPriorProbSpliceAcceptorSNS(variant, boundaries, variantData):
    '''
    Given a variant, boundaries (either PRIORS or ENIGMA), and list of dictionaries with variant data
    Determines reference and de novo acceptor scores for variant
      If variant in exon, also determines de novo donor scores and protein prior
    If variant causes a nonsense mutation, determines if splice rescue occurs
    Returns dicitionary containing scores for ref and de novo splice donor/acceptor
        score = "-" if score not applicable for variant
    Also contains other values:
        applicable prior, highest prior if variant has multiple priors
        applicable classe, highest predicted qualitative enigma class based on applicable prior
        ref prior, prior prob for reference splice sequence
        de novo donor and acceptor priors, prior prob for de novo splice sequence
        splice flag = 1, because variant in reference splice site
        de novo acc flag = 1 if variant is possible de novo acceptor variant
        de novo donor flag = 1 if variant is possible de novo donor variant
        spliceRescue = 1 if splice rescue possible for nonsense variant, 0 otherwise
        spliceFlag = 1 if splice rescue is possible so variant can be flagged for further analysis, 0 otherwise
        frameshift = 1 if nonsense variant causes a frameshift mutation also, 0 otherwise
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
            deNovoDonorInfo = {"refMaxEntScanScore": "-",
                               "refZScore": "-",
                               "altMaxEntScanScore": "-",
                               "altZScore": "-",
                               "deNovoDonorFlag": 0,
                               "priorProb": "N/A"}

        spliceRescue = 0
        spliceFlag = 0
        frameshift = 0
        # to check for nonsense variants in exonic portion of splice acceptor site
        if varInExon(variant) == True and getVarConsequences(variant) == "stop_gained":
            nonsenseData = getPriorProbSpliceRescueNonsenseSNS(variant, boundaries)
            applicablePrior = nonsenseData["priorProb"]
            spliceRescue = nonsenseData["spliceRescue"]
            spliceFlag = nonsenseData["spliceFlag"]
            frameshift = nonsenseData["frameshift"]
        
        return {"applicablePrior": applicablePrior,
                "applicableEnigmaClass": getEnigmaClass(applicablePrior),
                "proteinPrior": proteinPrior,
                "refDonorPrior": "-",
                "deNovoDonorPrior": deNovoDonorPrior,
                "refRefDonorMES": "-",
                "refRefDonorZ": "-",
                "altRefDonorMES": "-",
                "altRefDonorZ": "-",
                "refDeNovoDonorMES": deNovoDonorInfo["refMaxEntScanScore"],
                "refDeNovoDonorZ": deNovoDonorInfo["refZScore"],
                "altDeNovoDonorMES": deNovoDonorInfo["altMaxEntScanScore"],
                "altDeNovoDonorZ": deNovoDonorInfo["altZScore"],
                "deNovoDonorFlag": deNovoDonorInfo["deNovoDonorFlag"],
                "deNovoAccPrior": deNovoAccInfo["priorProb"],
                "refAccPrior": refPrior,
                "refRefAccMES": refSpliceInfo["refMaxEntScanScore"],
                "refRefAccZ": refSpliceInfo["refZScore"],
                "altRefAccMES": refSpliceInfo["altMaxEntScanScore"],
                "altRefAccZ": refSpliceInfo["altZScore"],
                "refDeNovoAccMES": deNovoAccInfo["refMaxEntScanScore"],
                "refDeNovoAccZ": deNovoAccInfo["refZScore"],
                "altDeNovoAccMES": deNovoAccInfo["altMaxEntScanScore"],
                "altDeNovoAccZ": deNovoAccInfo["altZScore"],
                "deNovoAccFlag": deNovoAccInfo["deNovoAccFlag"],
                "spliceSite": refSpliceInfo["spliceSite"],
                "spliceRescue": spliceRescue,
                "spliceFlag": spliceFlag,
                "frameshift": frameshift}
    
def getPriorProbProteinSNS(variant, variantData):
    '''
    Given a variant and a list of dictionaries containing variant data,
    Returns a dictionary containing:
      the variant's protein prior probability and enigma class for that prior
    '''
    if getVarType(variant) == "substitution":
        varHGVS = variant["HGVS_cDNA"]
        varGene = variant["Gene_Symbol"]

        for var in variantData:
            if var['gene'] == varGene and var['nthgvs'] == varHGVS:
                proteinPrior = float(var["protein_prior"])
                return {"priorProb": proteinPrior,
                        "enigmaClass": getEnigmaClass(proteinPrior)}

def getPriorProbProteinSNS(variant, variantData):
    '''
    Given a variant and a file containing protein prior probabilities,
    Returns a dictionary containing:
      the variant's protein prior probability and enigma class for that prior
    '''
    if getVarType(variant) == "substitution":
        varHGVS = variant["HGVS_cDNA"]
        varGene = variant["Gene_Symbol"]
        
        for var in variantData:
            if var['gene'] == varGene and var['nthgvs'] == varHGVS:
                proteinPrior = float(var["protein_prior"])
                return {"priorProb": proteinPrior,
                        "enigmaClass": getEnigmaClass(proteinPrior)}

def getPriorProbInGreyZoneSNS(variant, boundaries, variantData):
    '''
    Given a variant and a variant file, return prior prob and enigma class
    '''
    if getVarType(variant) == "substitution" and getVarLocation(variant, boundaries) == "grey_zone_variant":
        proteinData = getPriorProbProteinSNS(variant, variantData)
        return {"priorProb": proteinData["priorProb"],
                "enigmaClass": proteinData["enigmaClass"]}
                
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
                  "refRefDonorMES", "refRefDonorZ", "altRefDonorMES", "altRefDonorZ", "refDeNovoDonorMES", "refDeNovoDonorZ", "altDeNovoDonorMES",
                  "altDeNovoDonorZ", "deNovoDonorFlag", "refAccPrior", "deNovoAccPrior", "refRefAccMES", "refRefAccZ", "altRefAccMES",
                  "altRefAccZ", "refDeNovoAccMES", "refDeNovoAccZ", "altDeNovoAccMES", "altDeNovoAccZ", "deNovoAccFlag", "spliceSite",
                  "spliceRescue", "spliceFlag", "frameshift"]
    for header in newHeaders:
        fieldnames.append(header)
    outputData = csv.DictWriter(open(args.outputFile, "w"), delimiter="\t", fieldnames=fieldnames)
    outputData.writerow(dict((fn,fn) for fn in inputData.fieldnames))
    
    for variant in inputData:
        varDict = getVarDict(variant)
    # TO DO - create conditional to account for user selected boundaries
    # TO DO - create built_with_priors (copy of built)
    
if __name__ == "__main__":
    main()
