import json
import os
import re
import subprocess

from Bio.Seq import Seq
from pyfaidx import Fasta

from calc_priors.constants import BRCA1_CANONICAL, BRCA2_CANONICAL, BRCA_ZSCORES
from calc_priors import verify

import calcMaxEntScanMeanStd


def getExonBoundaries(variant):
    """
    Given a variant, returns the exon boundaries for the variant's transcript in a dictionary with format:
    key = exon number, value = dictionary with exon start and exon end for specific exon
    Uses function implemented in calcMaxEntScanMeanStd to get data for variant's transcript
    """
    varTranscript = variant["Reference_Sequence"]
    transcriptData = verify.getTranscriptData(varTranscript)
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


def getVarChrom(variant):
    """Given a variant, returns the chromosome based on the variant's gene_symbol"""
    varGene = variant["Gene_Symbol"]

    if varGene == "BRCA1":
        return "chr17"
    elif varGene == "BRCA2":
        return "chr13"
    else:
        return ""


def getVarConsequences(variant):
    """
    Given a variant, uses Ensembl VEP API to get variant consequences
    (e.g. intron variant, frameshift variant, missense variant)
    using variant chromosome, Hg38 start, Hg38 end, and alternate allele as input for API
    returns a list of strings detailing consequences of variant
    """

    # varStrand always 1 because all alternate alleles and positions refer to the plus strand
    varStrand = 1
    varAlt = variant["Alt"]

    assert variant["Chr"] in ["13", "17"]
    for base in varAlt:
        # API only works for alt alleles that are composed of the 4 canonical bases
        assert base in ["A", "C", "G", "T"]

    query = "%s:%s-%s:%s/%s" % (variant["Chr"], variant["Hg38_Start"],
                                variant["Hg38_End"], varStrand, varAlt)
    # Query local vep using query minus '?' character
    cmd = ["vep", "--cache", "--dir_cache", "/references/vep/",
           "--no_stats", "--offline", "--fasta", "/references/hg38.fa",
           "--output_file", "STDOUT", "--json", "--input_data", query]
    vep = json.loads(subprocess.check_output(cmd))

    # Should only be one BRCA1 canonical transcript in the list
    assert len([gene["consequence_terms"] for gene in vep["transcript_consequences"]
                if gene["transcript_id"] == BRCA1_CANONICAL
                or gene["transcript_id"] == BRCA2_CANONICAL]) == 1

    for gene in vep["transcript_consequences"]:
        if gene["transcript_id"] == BRCA1_CANONICAL or gene["transcript_id"] == BRCA2_CANONICAL:
            return gene["consequence_terms"]


def getVarType(variant):
    """
    Returns a string describing type of variant
    -substitution, deletion, insertion, delins, other
    depending on variant reference and alternate alleles
    """
    varRef = variant["Ref"]
    varAlt = variant["Alt"]
    acceptableRefSeq = verify.checkSequence(varRef)
    acceptableAltSeq = verify.checkSequence(varAlt)

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


def getRefSpliceDonorBoundaries(variant, intronicLength, exonicLength):
    """
    Given a variant, intronicLength and exonicLength returns the splice donor boundaries
    intronicLength = number of bp in intron that will be considered as part of splice donor region
    exonicLength = number of bp in exon that will be considered as part of splice donor region
    splice region is the last exonicLength bp in the exon and first intronicLength bp in the intron
    for the variant's transcript in a dictionary with the format:
    key = exon number, value = dictionary with donor start and donor end for exon
    """
    varExons = getExonBoundaries(variant)
    donorExons = varExons.copy()
    if variant["Gene_Symbol"] == "BRCA1":
        del donorExons["exon24"]
    elif variant["Gene_Symbol"] == "BRCA2":
        del donorExons["exon27"]
    varStrand = verify.getVarStrand(variant)
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
    """
    Given a variant, intronicLength and exonicLength returns the splice acceptor boundaries
    intronicLength = number of bp in intron that will be considered as part of splice acceptor region
    exonicLength = number of bp in exon that will be considered as part of splice acceptor region
    splice rgion is the last intronicLength bp in the exon and first exonicLength bp in the exon
    for the variant's transcript in a dictionary with the format:
    key = exon number, value = a dictionary with acceptor start and acceptor end for exon
    """
    varExons = getExonBoundaries(variant)
    acceptorExons = varExons.copy()
    if variant["Gene_Symbol"] == "BRCA1" or variant["Gene_Symbol"] == "BRCA2":
        del acceptorExons["exon1"]
    varStrand = verify.getVarStrand(variant)
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


def getFastaSeq(chrom, rangeStart, rangeStop, plusStrandSeq=True):
    """
    Given chromosome (in format 'chr13'), region genomic start position, and
    region genomic end position:
    Returns a string containing the sequence inclusive of rangeStart and rangeStop
    If plusStrandSeq=True, returns plus strand sequence
    If plusStrandSeq=False, returns minus strand sequence
    """
    if rangeStart < rangeStop:
        regionStart = rangeStart
        regionEnd = rangeStop
    else:
        regionStart = rangeStop
        regionEnd = rangeStart

    # NOTE: pyfaidx is NOT thread safe. Would be better to have one
    # REMIND: Switch to one per thread/process
    hg38 = Fasta("/references/hg38.fa", sequence_always_upper=True)
    sequence = hg38[chrom][regionStart - 1:regionEnd]

    if plusStrandSeq:
        return sequence.seq
    else:
        return sequence.reverse.complement.seq


def getSeqLocDict(chrom, varStrand, rangeStart, rangeStop):
    """
    Given chromosome, strand, region genomic start position, and region genomic end position
    returns a dictionary containing the genomic position as the key and reference allele as the value
    For minus strand gene (rangeStart > rangeStop), for plus strand gene (rangeStart < rangeStop)
    Always returns plus strand sequence
    """
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
    """
    Given a variant and a dictionary containing a sequence with bases and their locations,
    returns the dictionary with the alternate allele in place of the reference allele
    at the variant's genomic position
    """
    varRef = variant["Ref"]
    varAlt = variant["Alt"]
    varGenPos = int(variant["Pos"])
    seqLocDictRef = seqLocDict[varGenPos]
    if varRef == seqLocDictRef:
        altSeqDict = seqLocDict.copy()
        altSeqDict[varGenPos] = varAlt
    return altSeqDict


def getAltSeq(altSeqDict, varStrand):
    """
    Given a dictionary containing an alternate sequence with bases and their locations
    and the strand that the alternate allele is on
    Returns a string of the sequence containing the alternate allele
    """
    sequence = ""
    # to ensure that items in dictionary are sorted numerically
    for key, value in sorted(altSeqDict.items()):
        sequence += altSeqDict[key]
    if varStrand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    return sequence


def getRefAltSeqs(variant, rangeStart, rangeStop):
    """
    Given a variant, rangeStart, and rangeStop:
    Returns a dicitonary with ref and alt seq for the specified variant and range
    """
    varChrom = getVarChrom(variant)
    varStrand = verify.getVarStrand(variant)
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
    """
    Given a reference sequence and alternate sequence for a substitution variant
    Determines the index at which the alt seq differs from the ref seq
      - using zero indexing
    Returns the index at which altSeq differs from refSeq
    For example if the refSeq is ACTG and altSeq is AGTG, returns 1
    """
    if len(refSeq) == len(altSeq):
        for index in xrange(len(refSeq)):
            if refSeq.lower()[index] != altSeq.lower()[index]:
                return index
            else:
                pass
    else:
        return "N/A"


def getZScore(maxEntScanScore, donor=False):
    """
    Given a MaxEntScanScore, returns the zscore
    If donor is True, uses splice donor mean and std
    If donor is False, uses splice acceptor mean and std
    """
    # stdMeanData = json.load(open(os.path.join(os.path.dirname(__file__), 'brca.zscore.json'), "r"))
    stdMeanData = BRCA_ZSCORES

    if not donor:
        std = stdMeanData["acceptors"]["std"]
        mean = stdMeanData["acceptors"]["mean"]
    else:
        std = stdMeanData["donors"]["std"]
        mean = stdMeanData["donors"]["mean"]

    zscore = (maxEntScanScore - mean) / std
    return zscore


def getRefAltScores(refSeq, altSeq, donor=False):
    """
    Given ref and alt sequences and if sequence is in a splice donor region or not (True/False)
    Returns a dictionary containing raw MaxEntScan scores and zscores for ref and alt sequences
    """
    if not donor:
        refMaxEntScanScore = calcMaxEntScanMeanStd.runMaxEntScan(refSeq, donor=False)
        refZScore = getZScore(refMaxEntScanScore, donor=False)
        altMaxEntScanScore = calcMaxEntScanMeanStd.runMaxEntScan(altSeq, donor=False)
        altZScore = getZScore(altMaxEntScanScore, donor=False)
    else:
        refMaxEntScanScore = calcMaxEntScanMeanStd.runMaxEntScan(refSeq, donor=True)
        refZScore = getZScore(refMaxEntScanScore, donor=True)
        altMaxEntScanScore = calcMaxEntScanMeanStd.runMaxEntScan(altSeq, donor=True)
        altZScore = getZScore(altMaxEntScanScore, donor=True)

    scoreDict = {"refScores": {"maxEntScanScore": refMaxEntScanScore,
                               "zScore": refZScore},
                 "altScores": {"maxEntScanScore": altMaxEntScanScore,
                               "zScore": altZScore}}
    return scoreDict


def getMaxEntScanScoresSlidingWindowSNS(variant, windowSize, donor=False):
    """
    Given a variant and window size determines window sequences and scores for a sliding window
      that is the size of windowSize
    If donor=True, calculates MaxEntScan scores for splice donors
    If donor=False, calculates MaxEntScan scores for splice acceptors
    Returns a dictionary containing:
        1. window sequences - ref and alt seq for each window (variant in positions 1-windowSize)
        2. window scores - ref and alt MaxEntScan scores and zscores for each window
        3. window alt MaxEntScan scores - only contains alt MaxEntScan scores for each window
    """
    varGenPos = int(variant["Pos"])
    varStrand = verify.getVarStrand(variant)
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


def getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, exonicPortionSize, intronicPortionSize,
                         donor=True):
    """
    Given a variant's:
    genetic postion, strand, sliding window position with max MES score AND
      whether that position is within exonic portion of highest scoring window, exonic portion size, and intronic portion size
    Returns the position where splicing occurs for a de novo splice donor or acceptor (depending on donor=True argument)
    """
    if varStrand == "+":
        if not inExonicPortion:
            if donor:
                newSplicePos = int(varGenPos) - (varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) + abs(varWindowPos - intronicPortionSize)
        else:
            if donor:
                newSplicePos = int(varGenPos) + abs(varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) - (varWindowPos - intronicPortionSize)
    else:
        if not inExonicPortion:
            if donor:
                newSplicePos = int(varGenPos) + (varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) - abs(varWindowPos - intronicPortionSize)
        else:
            if donor:
                newSplicePos = int(varGenPos) - abs(varWindowPos - exonicPortionSize)
            else:
                newSplicePos = int(varGenPos) + (varWindowPos - intronicPortionSize)
    return newSplicePos


def getEnigmaClass(priorProb):
    """
    Given a prior probability of pathogenecity, returns a predicted qualitative ENIGMA class
    """
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


