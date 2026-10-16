import hgvs.dataproviders.uta as uta
import hgvs.assemblymapper as am
import hgvs.parser

from calcMaxEntScanMeanStd import fetch_gene_coordinates
from calc_priors.constants import brca1CIDomains, brca2CIDomains, greyZones, BRCA1_RefSeq, \
    BRCA2_RefSeq

_CHROM_TO_NC = {
    'chr13': 'NC_000013.11',
    'chr17': 'NC_000017.11',
}
_hp = hgvs.parser.Parser()
_mapper = None  # initialised lazily on first call so UTA_DB_URL can be set at runtime


def _get_mapper():
    global _mapper
    if _mapper is None:
        _hdp = uta.connect()
        _mapper = am.AssemblyMapper(_hdp, assembly_name='GRCh38', alt_aln_method='splign',
                                    prevalidation_level=None)
    return _mapper

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
    if transcriptData is None:
        return True
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
    Given a reference sequence (e.g. "NM_007294.3" or "NM_007294.4"),
    Returns transcript data for that reference sequence.
    Version suffix is ignored so that newer RefSeq versions still resolve.
    """
    base = referenceSequence.split('.')[0] if referenceSequence else ''
    if base == BRCA1_RefSeq.split('.')[0]:
        return brca1TranscriptData
    elif base == BRCA2_RefSeq.split('.')[0]:
        return brca2TranscriptData

    print("Warning: No reference sequence available for " + str(referenceSequence))
    return None


def convertGenomicPosToTranscriptPos(genomicPos, chrom, transcript):
    """
    Given a genomic position, chrom (in format "chrN"), and transcript accession string:
    Returns a string of the transcript (cDNA) position at the given genomic position.
    Uses the biocommons hgvs library via UTA for coordinate mapping.
    """
    chrom_key = chrom if chrom.startswith('chr') else f'chr{chrom}'
    nc = _CHROM_TO_NC[chrom_key]
    g_var = _hp.parse_hgvs_variant(f'{nc}:g.{genomicPos}T>A')
    c_var = _get_mapper().g_to_c(g_var, transcript)
    return str(c_var.posedit.pos.start)


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
