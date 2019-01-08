import collections


# idea from https://stackoverflow.com/a/28452633/346905
from calc_priors.compute import getVarLocation
from calc_priors.extract import getVarType
from calc_priors.verify import getVarStrand


class ReadOnlyDict(collections.Mapping):
    def __init__(self, data):
        self._data = data

    def __getitem__(self, key):
        return self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)


# FIXME: make this read-only again and update addVarDataToRow to base its copy off this
BLANK_DICT = {
    "applicablePrior": "-",
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
     "lowMESFlag": "-"
}


def addVarDataToRow(varData, inputRow):
    """
    Given data about a particular variant and a row from input file,
    Returns row with appended data
    """
    for key in varData.keys():
        inputRow[key] = varData[key]
    return inputRow


def getVarDict(variant, boundaries):
    """
    Given input data, returns a dictionary containing information for each variant in input
    Dictionary key is variant HGVS_cDNA and value is a dictionary containing variant gene, variant chromosome,
    variant strand, variant genomic coordinate, variant type, and variant location
    """
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
