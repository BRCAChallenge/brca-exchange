import pyhgvs

from calcMaxEntScanMeanStd import fetch_gene_coordinates
from calc_priors.constants import brca1CIDomains, brca2CIDomains, greyZones, BRCA1_RefSeq, \
    BRCA2_RefSeq

# Fetch transcript data for BRCA1/BRCA2 RefSeq transcripts
brca1TranscriptData = fetch_gene_coordinates(BRCA1_RefSeq)
brca2TranscriptData = fetch_gene_coordinates(BRCA2_RefSeq)


def checkSequence(sequence):
    """Checks if a given sequence contains acceptable nucleotides returns True if sequence is comprised entirely of acceptable bases"""
    acceptableBases = ["A", "C", "T", "G", "N", "R", "Y"]
    if len(sequence) > 0:
        for base in sequence:
            if base not in acceptableBases:
                return False
        return True
    else:
        return False


def checkWithinBoundaries(varStrand, varGenPos, boundaryStart, boundaryEnd):
    """
    Checks whether a position (varGenPos) is within certain boundaries (boundaryStart and boundaryEnd)
    Dependent on the variant's transcript strand (varStrand)
    Function is inclusive of boundaryStart and boundaryEnd
    """
    if varStrand == "+":
        if varGenPos >= boundaryStart and varGenPos <= boundaryEnd:
            return True
    elif varStrand == "-":
        if varGenPos <= boundaryStart and varGenPos >= boundaryEnd:
            return True
    else:
        return False


def varOutsideBoundaries(variant):
    """Given a variant, determines if variant is outside transcript boundaries"""
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
    """
    Given a variant, if variant is inside transcript boundaries,
    determines if variant is in 3' or 5' UTR of transcript
    """
    varOutBounds = varOutsideBoundaries(variant)
    if not varOutBounds:
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


def varInGreyZone(variant):
    """
    Given a variant, determines if variant is in the grey zone
    Returns True if variant is in the grey zone
    Returns False if variant is NOT in the grey zone
    """
    varGenPos = int(variant["Pos"])
    varGene = variant["Gene_Symbol"]
    varStrand = getVarStrand(variant)
    if varGene == "BRCA2":
        greyZoneStart = greyZones[varGene]["greyZoneStart"]
        greyZoneEnd = greyZones[varGene]["greyZoneEnd"]
        withinBoundaries = checkWithinBoundaries(varStrand, varGenPos, greyZoneStart, greyZoneEnd)
        if withinBoundaries:
            return True
    return False


def varAfterGreyZone(variant):
    """
    Given a variant, determines if variant is after the gene grey zone
    Returns True if variant is after the grey zone
    """
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


def isCIDomainInRegion(regionStart, regionEnd, boundaries, gene):
    """
    Given a region of interest, boundaries (either enigma or priors) and gene of interest
    Determines if there is an overlap between the region of interest and a CI domain
    For minus strand gene (BRCA1) regionStart > regionEnd
    For plus strand gene (BRCA2) regionStart < regionEnd
    Returns True if there is an overlap, False otherwise
    """
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


def compareRefAltExonLengths(refLength, altLength):
    """
    Compares ref and alt exon lengths
    If both exon lengths % 3 are equal, then both ref and alt have the same reading frame
    Returns true if both ref and alt exon have the same reading frame, false otherwise
    """
    if refLength % 3 == altLength % 3:
        return True
    else:
        return False


# FIXME: getVarStrand(),  convertGenomicPosToTranscriptPos(), and formatSplicePosition() aren't predicates
# they should be moved to another module, but currently the most logical place (extract) creates a circular dependency
# between extract and verify, since those methods are used in this module and verify is used in extract

def getVarStrand(variant):
    """Given a variant, returns the coding strand based on the variant's gene_symbol"""
    varGene = variant["Gene_Symbol"]

    if varGene == "BRCA1":
        return '-'
    elif varGene == "BRCA2":
        return '+'
    else:
        return ""


def getTranscriptData(referenceSequence):
    """
    Given a reference sequence (e.g. "NM_007294.3"),
    Returns transcript data for that reference sequencee
    """
    if referenceSequence == BRCA1_RefSeq:
        return brca1TranscriptData
    elif referenceSequence == BRCA2_RefSeq:
        return brca2TranscriptData


def convertGenomicPosToTranscriptPos(genomicPos, chrom, genome, transcript):
    """
    Given a genomic position, chrom (in format "chrN"), genome (SequenceFileDB for genome),
      and transcript (pyhgvs transcript object):
    Returns a string of the transcript position at the given genomic position
    """
    # use "T" and "A" for ref and alt because transcript position is not dependent on these values
    # converts genomic position to transcript position
    hgvs_name = str(pyhgvs.format_hgvs_name(chrom, genomicPos, "T", "A", genome, transcript))
    # parses out transcript position from full hgvs_name
    transcriptPos = str(pyhgvs.HGVSName(hgvs_name).cdna_start)
    return transcriptPos


def formatSplicePosition(position, transcript=False):
    """
    Given a position and transcript argument, returns a formatted splice position
    If transcript = True, returns transcript formatted position "c.N"
    If transcript = False, returns genomic formatted position "g.N"
    """
    if transcript:
        return "c." + str(position)
    else:
        return "g." + str(position)