import unittest
import mock
import calcVarPriors
import calcMaxEntScanMeanStd
from calcVarPriorsMockedResponses import brca1Exons, brca2Exons 
from calcVarPriorsMockedResponses import brca1RefSpliceDonorBounds, brca2RefSpliceDonorBounds 
from calcVarPriorsMockedResponses import brca1RefSpliceAcceptorBounds, brca2RefSpliceAcceptorBounds

# dictionary containing possible strands for variants
strand = {"minus": "-",
          "plus": "+"}

# dictionary containing possible chromosomes for variants
chromosomes = {"13": "chr13",
               "17": "chr17"}

# dictionary containing possible variant types
varTypes = {"sub": "substitution",
            "ins": "insertion",
            "del": "deletion",
            "delins": "delins"}

# dictionary containing the number of exons for particular genes
numExons = {"BRCA1": 23,
            "BRCA2": 27}

# example donor exon boundaries for BRCA1 and BRCA2
exonDonorBoundsBRCA1 = {"exon16": {"donorStart": 43070930,
                                   "donorEnd": 43070922}}
exonDonorBoundsBRCA2 = {"exon15": {"donorStart": 32356607,
                                   "donorEnd": 32356615}}

# example acceptor exon boundaries for BRCA1 and BRCA2
exonAcceptorBoundsBRCA1 = {"exon21": {"acceptorStart": 43051137,
                                      "acceptorEnd": 43051115}}
exonAcceptorBoundsBRCA2 = {"exon20": {"acceptorStart": 32370936,
                                      "acceptorEnd": 32370958}}

# dictionary containing possible variant locations
variantLocations = {"outBounds": "outside_transcript_boundaries_variant",
                    "inCI": "CI_domain_variant",
                    "inCISpliceDonor": "CI_splice_donor_variant",
                    "inCISpliceAcceptor": "CI_splice_acceptor_variant",
                    "inSpliceDonor": "splice_donor_variant",
                    "inSpliceAcceptor": "splice_acceptor_variant",
                    "inGreyZone": "grey_zone_variant",
                    "afterGreyZone": "after_grey_zone_variant",
                    "inExon": "exon_variant",
                    "inUTR": "UTR_variant",
                    "inIntron": "intron_variant"}

# transcript data to mock response from fetch_gene_coordinates for BRCA1 and BRCA2
transcriptDataBRCA1 = {'bin': '114',
                       'exonEnds': '43045802,43047703,43049194,43051117,43057135,43063373,43063951,43067695,43071238,43074521,43076614,43082575,43091032,43094860,43095922,43097289,43099880,43104261,43104956,43106533,43115779,43124115,43125483,',
                       'exonFrames': '1,0,1,0,0,1,1,0,1,2,1,0,1,1,2,1,0,1,2,2,2,0,-1,',
                       'name': 'NM_007294.3',
                       'txStart': 43044294,
                       'exonCount': 23,
                       'cdsEndStat': 'cmpl',
                       'cdsEnd': 43124096,
                       'score': 0,
                       'name2': 'BRCA1',
                       'strand': '-',
                       'cdsStart': 43045677,
                       'cdsStartStat': 'cmpl',
                       'chrom': 'chr17',
                       'txEnd': 43125483,
                       'exonStarts': '43044294,43047642,43049120,43051062,43057051,43063332,43063873,43067607,43070927,43074330,43076487,43082403,43090943,43091434,43095845,43097243,43099774,43104121,43104867,43106455,43115725,43124016,43125270,'}

transcriptDataBRCA2 = {'bin': '103',
                       'exonEnds': '32315667,32316527,32319325,32325184,32326150,32326282,32326613,32329492,32331030,32333387,32341196,32344653,32346896,32355288,32356609,32357929,32362693,32363533,32370557,32371100,32376791,32379515,32379913,32380145,32394933,32397044,32399672,',
                       'exonFrames': '-1,0,1,1,2,1,0,1,0,1,1,1,1,2,1,0,2,2,0,0,1,0,1,0,1,0,0,',
                       'name': 'NM_000059.3',
                       'txStart': 32315479,
                       'exonCount': 27,
                       'cdsEndStat': 'cmpl',
                       'cdsEnd': 32398770,
                       'score': 0,
                       'name2': 'BRCA2',
                       'strand': '+',
                       'cdsStart': 32316460,
                       'cdsStartStat': 'cmpl',
                       'chrom': 'chr13',
                       'txEnd': 32399672,
                       'exonStarts': '32315479,32316421,32319076,32325075,32326100,32326241,32326498,32329442,32330918,32332271,32336264,32344557,32346826,32354860,32356427,32357741,32362522,32363178,32370401,32370955,32376669,32379316,32379749,32380006,32394688,32396897,32398161,'}

# sample sequence data for BRCA1 and BRCA2
brca1Seq = "GATCTGGAAGAAGAGAGGAAGAG"
brca2Seq = "TGTGTAACACATTATTACAGTGG"

# MaxEntScan score mean and std for donors and acceptors
meanStdDict =  {"donors": {"std": 2.3289956850167082,
                           "mean": 7.9380909090909073},
                "acceptors": {"std": 2.4336623152078452,
                              "mean": 7.984909090909091}}

# possible predicted qualitative ENIGMA classes
enigmaClasses = {"class1": "class_1",
                 "class2": "class_2",
                 "class3": "class_3",
                 "class4": "class_4",
                 "class5": "class_5"}

# possible prior probability of pathogenecity values
priorProbs = {"low": 0.04,
              "moderate": 0.34,
              "high": 0.97,
              "NA": "N/A"}


class test_calcVarPriors(unittest.TestCase):

    def setUp(self):

        self.variant = {"Chr": "13",
                        "Pos": "32314943",
                        "Ref": "A",
                        "Alt": "G",
                        "Gene_Symbol": "BRCA2",
                        "Reference_Sequence": "NM_000059.3",
                        "pyhgvs_cDNA": "NM_000059.3:c.-764A>G"}
                  
    def test_checkSequence(self):
        '''Tests that checkSequence function categorized acceptable sequences correctly'''
        # sequence with unacceptable letters
        self.variant["Ref"] = "ATGSFHG"
        self.variant["Alt"] = "AGTHA"
        acceptableRefSeq = calcVarPriors.checkSequence(self.variant["Ref"])
        acceptableAltSeq = calcVarPriors.checkSequence(self.variant["Alt"])
        self.assertFalse(acceptableRefSeq)
        self.assertFalse(acceptableAltSeq)

        # sequence with numbers
        self.variant["Ref"] = "3452345"
        self.variant["Alt"] = "3456324"
        acceptableRefSeq = calcVarPriors.checkSequence(self.variant["Ref"])
        acceptableAltSeq = calcVarPriors.checkSequence(self.variant["Alt"])
        self.assertFalse(acceptableRefSeq)
        self.assertFalse(acceptableAltSeq)

        # blank sequence
        self.variant["Ref"] = ""
        self.variant["Alt"] = ""
        acceptableRefSeq = calcVarPriors.checkSequence(self.variant["Ref"])
        acceptableAltSeq = calcVarPriors.checkSequence(self.variant["Alt"])
        self.assertFalse(acceptableRefSeq)
        self.assertFalse(acceptableAltSeq)

        # sequence with only ATCG
        self.variant["Ref"] = "ATGACG"
        self.variant["Alt"] = "AGTAATA"
        acceptableRefSeq = calcVarPriors.checkSequence(self.variant["Ref"])
        acceptableAltSeq = calcVarPriors.checkSequence(self.variant["Alt"])
        self.assertTrue(acceptableRefSeq)
        self.assertTrue(acceptableAltSeq)

        # sequence containing all possible acceptable bases
        self.variant["Ref"] = "ATGRACYGN"
        self.variant["Alt"] = "YAGRTNAATA"
        acceptableRefSeq = calcVarPriors.checkSequence(self.variant["Ref"])
        acceptableAltSeq = calcVarPriors.checkSequence(self.variant["Alt"])
        self.assertTrue(acceptableRefSeq)
        self.assertTrue(acceptableAltSeq)

    def test_getVarStrand(self):
        '''Tests that variant strand is set correctly based on variant's gene_symbol'''
        self.variant["Gene_Symbol"] = "BRCA1"
        varStrand = calcVarPriors.getVarStrand(self.variant)
        self.assertEquals(varStrand, strand["minus"])

        self.variant["Gene_Symbol"] = "BRCA2"
        varStrand = calcVarPriors.getVarStrand(self.variant)
        self.assertEquals(varStrand, strand["plus"])

    def test_getVarChrom(self):
        '''Tests taht variant chromosome is set correctly based on variant's gene_symbol'''
        self.variant["Gene_Symbol"] = "BRCA1"
        varChrom = calcVarPriors.getVarChrom(self.variant)
        self.assertEquals(varChrom, chromosomes["17"])

        self.variant["Gene_Symbol"] = "BRCA2"
        varChrom = calcVarPriors.getVarChrom(self.variant)
        self.assertEquals(varChrom, chromosomes["13"])

    @mock.patch('calcVarPriors.checkSequence', return_value = True)
    def test_getVarType(self, checkSequence):
        '''
        Tests that variant type is set correctly to substitution, deletion, insertion, or delins based on variant "Ref" and "Alt" values
        '''
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, varTypes["sub"])

        self.variant["Ref"] = "A"
        self.variant["Alt"] = "AAA"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, varTypes["ins"])

        self.variant["Ref"] = "AGT"
        self.variant["Alt"] = "A"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, varTypes["del"])

        self.variant["Ref"] = "AG"
        self.variant["Alt"] = "AGTA"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, varTypes["delins"])

        self.variant["Ref"] = "AGTA"
        self.variant["Alt"] = "AG"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, varTypes["delins"])

        self.variant["Ref"] = "AG"
        self.variant["Alt"] = "GT"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, varTypes["delins"])

    def test_getVarConsequences(self):
        '''
        Tests that:
        1. Variants with non-BRCA1/BRCA2 chromosomes are skipped
        2. Variants with Alt alleles that are not one of the 4 canonical bases are skipped
        '''

        self.variant["Chr"] = ""
        varCons = calcVarPriors.getVarConsequences(self.variant)
        self.assertEquals(varCons, "unable_to_determine")

        self.variant["Chr"] = "41160094"
        varCons = calcVarPriors.getVarConsequences(self.variant)
        self.assertEquals(varCons, "unable_to_determine")

        self.variant["Chr"] = "chr17:g.43008077:TAGG"
        varCons = calcVarPriors.getVarConsequences(self.variant)
        self.assertEquals(varCons, "unable_to_determine")

        self.variant["Chr"] = "13"
        self.variant["Hg38_Start"] = "32339320"
        self.variant["Hg38_End"] = "32339320"
        self.variant["Alt"] = "R"
        varCons = calcVarPriors.getVarConsequences(self.variant)
        self.assertEquals(varCons, "unable_to_determine")

        self.variant["Alt"] = "-"
        varCons = calcVarPriors.getVarConsequences(self.variant)
        self.assertEquals(varCons, "unable_to_determine")

        self.variant["Alt"] = "38413620"
        varCons = calcVarPriors.getVarConsequences(self.variant)
        self.assertEquals(varCons, "unable_to_determine")

    def test_checkWithinBoundaries(self):
        '''
        Tests that positions are correctly identified as in/not in boundaries and that boundaries are inclusive
        '''
        varStrand = strand["plus"]
        boundaryStart = 32357742
        boundaryEnd = 32357780

        # check that last position before boundary start is NOT included
        position = 32357741
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertFalse(withinBoundaries)

        # check that first position after boundary end is NOT included
        position = 32357781
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertFalse(withinBoundaries)

        # check that boundaryStart is included
        position = 32357742
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertTrue(withinBoundaries)

        # check that boundaryEnd is included
        position = 32357780
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertTrue(withinBoundaries)

        # check that position within boundaries is included
        position = 32357758
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertTrue(withinBoundaries)

        varStrand = strand["minus"]
        boundaryStart = 43067695
        boundaryEnd = 43067649

        # check that last position before boundary start is NOT included
        position = 43067696
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertFalse(withinBoundaries)

        # check that first position after boundary end is NOT included
        position = 43067648
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertFalse(withinBoundaries)

        # check that boundaryStart is included
        position = 43067695
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertTrue(withinBoundaries)
        
        # check that boundaryEnd is included
        position = 43067649
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertTrue(withinBoundaries)

        # check that position within boundaries is included
        position = 43067669
        withinBoundaries = calcVarPriors.checkWithinBoundaries(varStrand, position, boundaryStart, boundaryEnd)
        self.assertTrue(withinBoundaries)        

    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_varOutsideBoundariesBRCA1(self, getVarStrand):
        '''Tests that variant outside/inside transcript boundaries are correctly identified for BRCA1'''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"

        # checks for BRCA1 variant outside transcript boundaries
        self.variant["Pos"] = "43044274"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertTrue(varOutBounds)

        # checks for BRCA1 variant inside transcript boundaries
        self.variant["Pos"] = "43070957"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertFalse(varOutBounds)

    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_varOutsideBoundariesBRCA2(self, getVarStrand):
        '''Tests that variant outside/inside transcript boundaries are correctly identified for BRCA2'''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"

        # checks for BRCA2 variant outside transcript boundaries
        self.variant["Pos"] = "32315477"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertTrue(varOutBounds)

        # checks for BRCA2 variant inside transcript boundaries
        self.variant["Pos"] = "32326500"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertFalse(varOutBounds)
        
    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_varInUTRBRCA1(self, varOutsideBoundaries, getVarStrand):
        '''Tests that variants in 5' and 3' UTR are correctly identified for BRCA1'''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        
        # checks for BRCA1 variant in 5' UTR
        self.variant["Pos"] = "43124110"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertTrue(varInUTR)
        
        # checks for BRCA1 variant in 3' UTR
        self.variant["Pos"] = "43045668"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertTrue(varInUTR)

        # checks for BRCA1 variant in exon
        self.variant["Pos"] = "43049184"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertFalse(varInUTR)

        # checks for BRCA1 variant in intron
        self.variant["Pos"] = "43049213"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertFalse(varInUTR)

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_varInUTRBRCA2(self, varOutsideBoundaries, getVarStrand):
        '''Tests that variants in 5' and 3' UTR are correctly identified for BRCA2'''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        
        # checks for BRCA2 variant in 5' UTR
        self.variant["Pos"] = "32316434"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertTrue(varInUTR)
        
        # checks for BRCA2 variant in 3' UTR
        self.variant["Pos"] = "32398781"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertTrue(varInUTR)

        # checks for BRCA2 variant in exon
        self.variant["Pos"] = "32394719"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertFalse(varInUTR)

        # checks for BRCA2 variant in intron
        self.variant["Pos"] = "32396875"
        varInUTR = calcVarPriors.varInUTR(self.variant)
        self.assertFalse(varInUTR)

    def test_getExonBoundariesBRCA1(self):
        '''
        Tests that:
        1. Exon boundaries are set correctly for gene on minus strand (BRCA1)
            - checks that exon does not start before it ends
            - next exon does not start before current exon ends
        2. length of varExons matches number of exons for BRCA1
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        varExons = calcVarPriors.getExonBoundaries(self.variant)
        # check correct number of exons
        self.assertEquals(len(varExons), numExons["BRCA1"])
        # because '-' strand gene, checks that exonStart > exonEnd
        for exon in varExons.keys():
            exonBounds = varExons[exon]
            self.assertGreater(exonBounds["exonStart"], exonBounds["exonEnd"])
            currentExonNum = int(exon[4:])
            nextExonNum = str(currentExonNum + 1)
            nextExonKey = "exon" + nextExonNum
            if nextExonNum <= numExons["BRCA1"]:
                nextExonStart = varExons[nextExonKey]["exonStart"]
                # checks that next exon does not start before current exon ends
                self.assertGreater(exonBounds["exonEnd"], nextExonStart)

    def test_getExonBoundariesBRCA2(self):
        '''
        Tests that:
        1. Exon boundaries are set correctly for gene on plus strand (BRCA2)
            - checks that exon does not start before it ends
            - next exon does not start before current exon ends
        2. length of varExons matches number of exons for BRCA2
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        varExons = calcVarPriors.getExonBoundaries(self.variant)
        # check correct number of exons
        self.assertEquals(len(varExons), numExons["BRCA2"])
        # because '+' strand gene, checks that exonEnd > exonStart
        for exon in varExons.keys():
            exonBounds = varExons[exon]
            self.assertGreater(exonBounds["exonEnd"], exonBounds["exonStart"])
            currentExonNum = int(exon[4:])
            nextExonNum = str(currentExonNum + 1)
            nextExonKey = "exon" + nextExonNum
            if nextExonNum <= numExons["BRCA2"]:
                nextExonStart = varExons[nextExonKey]["exonStart"]
                # checks that next exon does not start before current exon ends
                self.assertGreater(nextExonStart, exonBounds["exonEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_getRefSpliceDonorBoundariesBRCA1(self, getExonBoundaries, getVarStrand):
        '''
        Tests that splice donor boundaries are set correctly for reference transcript (NM_000059.3) and strand (-)
        Uses example exon boundaries for BRCA1 set in setUp function
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        spliceDonorBounds = calcVarPriors.getRefSpliceDonorBoundaries(self.variant)
        # checks that region after last exon is not considered a splice donor region
        self.assertNotIn("exon24", spliceDonorBounds)
        # to find exon specified in setUp function
        exon = exonDonorBoundsBRCA1.keys()[0]
        self.assertEquals(exonDonorBoundsBRCA1[exon]["donorStart"],
                          spliceDonorBounds[exon]["donorStart"])
        self.assertEquals(exonDonorBoundsBRCA1[exon]["donorEnd"],
                          spliceDonorBounds[exon]["donorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_getRefSpliceDonorBoundariesBRCA2(self, getExonBoundaries, getVarStrand):
        '''
        Tests that splice donor boundaries are set correctly for reference transcript (NM_000059.3) and strand (+)
        Uses example exon boundaries for BRCA2 set in setUp function
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        spliceDonorBounds = calcVarPriors.getRefSpliceDonorBoundaries(self.variant)
        # checks that region after last exon is not considered a splice donor region
        self.assertNotIn("exon27", spliceDonorBounds)
        # to find exon specified in setUp function
        exon = exonDonorBoundsBRCA2.keys()[0]
        self.assertEquals(exonDonorBoundsBRCA2[exon]["donorStart"],
                          spliceDonorBounds[exon]["donorStart"])
        self.assertEquals(exonDonorBoundsBRCA2[exon]["donorEnd"],
                          spliceDonorBounds[exon]["donorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_getRefSpliceAcceptorBoundariesBRCA1(self, getExonBoundaries, getVarStrand):
        '''
        Tests that splice acceptor boundaries are set correctly for reference transcript (NM_007294.3) and strand (-)
        Uses example exon boundaries for BRCA1 in setUp function
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        spliceAcceptorBounds = calcVarPriors.getRefSpliceAcceptorBoundaries(self.variant)
        # checks that region before first exon is not considered a splice acceptor region
        self.assertNotIn("exon1", spliceAcceptorBounds)
        # to find exon specified in setUp function
        exon = exonAcceptorBoundsBRCA1.keys()[0]
        self.assertEquals(exonAcceptorBoundsBRCA1[exon]["acceptorStart"],
                          spliceAcceptorBounds[exon]["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA1[exon]["acceptorEnd"],
                          spliceAcceptorBounds[exon]["acceptorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_getRefSpliceAcceptorBoundariesBRCA2(self, getExonBoundaries, getVarStrand):
        '''
        Tests that splice acceptor boundaries are set correctly for reference transcript (NM_000059.3) and strand (+)
        Uses example exon boundaries for BRCA2 in setUp function
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        spliceAcceptorBounds = calcVarPriors.getRefSpliceAcceptorBoundaries(self.variant)
        # checks that region before first exon is not considered a splice acceptor region
        self.assertNotIn("exon1", spliceAcceptorBounds)
        # to find exon specified in setUp function
        exon = exonAcceptorBoundsBRCA2.keys()[0]
        self.assertEquals(exonAcceptorBoundsBRCA2[exon]["acceptorStart"],
                          spliceAcceptorBounds[exon]["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA2[exon]["acceptorEnd"],
                          spliceAcceptorBounds[exon]["acceptorEnd"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_varInExonBRCA1(self, varOutsideBoundaries, getExonBoundaries, getVarStrand):
        '''Tests that variant is correctly identified as inside or outside an exon for BRCA1'''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"

        # checks BRCA1 5' exon boundary (last base in intron)
        self.variant["Pos"] = "43104957"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)

        # checks BRCA1 5' exon boundary (first base in exon)
        self.variant["Pos"] = "43104956"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # checks BRCA1 3' exon boundary (last base in exon)
        self.variant["Pos"] = "43067608"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # checks BRCA1 3' exon boundary (first base in intron)
        self.variant["Pos"] = "43067607"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)
        
        # checks BRCA1 variant inside an exon
        self.variant["Pos"] = "43049176"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # checks BRCA1 variant outside an exon
        self.variant["Pos"] = "43045827"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_varInExonBRCA2(self, varOutsideBoundaries, getExonBoundaries, getVarStrand):
        '''Tests that variant is correctly identified as inside or outside an exon for BRCA2'''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"

        # checks BRCA2 5' exon boundary (last base in intron)
        self.variant["Pos"] = "32357741"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)

        # checks BRCA2 5' exon boundary (first base in exon)
        self.variant["Pos"] = "32357742"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # checks BRCA2 3' exon boundary (last base in exon)
        self.variant["Pos"] = "32370557"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # checks BRCA2 3' exon boundary (first base in intron)
        self.variant["Pos"] = "32370558"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)
        
        # checks BRCA2 variant inside an exon
        self.variant["Pos"] = "32398201"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # cehcks BRCA2 variant outside an exon
        self.variant["Pos"] = "32396873"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)

    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca1RefSpliceDonorBounds)
    def test_varInSpliceRegionDonorBRCA1(self, getRefDonorBoundaries):
        '''
        Tests that:
        1. Variant is correctly identified as in or NOT in a splice donor region 
           for multiple positions across multiple exons in BRCA1
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        donor = True

        #checks that 7th base in intron is NOT counted as in splice donor
        self.variant["Pos"] = "43063326"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceDonor)
        
        # checks that first base in intron counted as in splice donor
        self.variant["Pos"] =  "43097243"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in intron counted as in splice donor
        self.variant["Pos"] = "43095842"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that last base in intron counted as in splice donor
        self.variant["Pos"] = "43091429"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that first base in exon counted as in splice donor
        self.variant["Pos"] = "43090946"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in exon counted as in splice donor
        self.variant["Pos"] = "43082405"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that last base in exon counted as in splice donor
        self.variant["Pos"] = "43076488"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that 4th to last base in exon is NOT counted as in splice donor
        self.variant["Pos"] = "43057055"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceDonor)

        # checks that region after exon 24 is  NOT counted as in splice donor
        self.variant["Pos"] = "43044294"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceDonor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca1RefSpliceDonorBounds)
    def test_getVarSpliceRegionBoundsDonorBRCA1(self, varInSpliceRegion, getRefSpliceDonorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position)
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        donor = True
    
        # checks that variant in exon 16 splice donor region boundaries are returned correctly
        self.variant["Pos"] = "43070925"
        spliceDonorRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor)
        self.assertEquals(exonDonorBoundsBRCA1["exon16"]["donorStart"], spliceDonorRegion["donorStart"])
        self.assertEquals(exonDonorBoundsBRCA1["exon16"]["donorEnd"], spliceDonorRegion["donorEnd"])

    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    def test_varInSpliceRegionDonorBRCA2(self, getRefSpliceDonorBoundaries):
        '''
        Tests that:
        1. Variant is correctly identified as in or NOT in a splice donor region 
           for multiple positions across multiple exons in BRCA2
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = True

        # checks that 7th base in intron is NOT counted as in splice donor
        self.variant["Pos"] = "32363540"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceDonor)

        # checks that first base in intron counted as in splice donor
        self.variant["Pos"] = "32329493"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in intron counted as in splice donor
        self.variant["Pos"] = "32331033"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that last base in intron counted as in splice donor
        self.variant["Pos"] = "32333393"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that first base in exon counted as in splice donor
        self.variant["Pos"] = "32341194"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in exon counted as in splice donor
        self.variant["Pos"] = "32344652"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that last base in exon counted as in splice donor
        self.variant["Pos"] = "32346896"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceDonor)

        # checks that 4th to last base in exon is NOT counted as in splice donor
        self.variant["Pos"] = "32370554"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceDonor)

        # checks that region after  exon 27 is NOT counted as in splice donor
        self.variant["Pos"] = "32399672"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceDonor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    def test_getVarSpliceRegionBoundsDonorBRCA2(self, varInSplcieRegion, getRefSpliceDonorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position) 
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = True

        # checks that variant in exon 16 splice donor region boundaries are returned correctly
        self.variant["Pos"] = "32356608"
        spliceDonorRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor)
        self.assertEquals(exonDonorBoundsBRCA2["exon15"]["donorStart"], spliceDonorRegion["donorStart"])
        self.assertEquals(exonDonorBoundsBRCA2["exon15"]["donorEnd"], spliceDonorRegion["donorEnd"])

    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_varInSpliceRegionAcceptorBRCA1(self, getRefSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Variant is correctly identified as in or NOT in a splice acceptor region
           uses multiple positions across multiple exons for BRCA1
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        donor = False

        # checks that -21st base in intron is NOT counted as in splice acceptor
        self.variant["Pos"] = "43067716"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceAcceptor)
        
        # checks that first base in intron counted as in splice acceptor
        self.variant["Pos"] = "43124135"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in intron counted as in splice acceptor
        self.variant["Pos"] = "43115787"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in intron counted as in splice acceptor
        self.variant["Pos"] = "43106534"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that first base in exon counted as in splice acceptor
        self.variant["Pos"] = "43104956"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in exon counted as in splice acceptor
        self.variant["Pos"] = "43104260"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in exon counted as in splice acceptor
        self.variant["Pos"] = "43099878"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that 4th base in exon is NOT counted as in splice acceptor
        self.variant["Pos"] = "43063948"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceAcceptor)

        # checks that region before exon 1 is NOT counted as in splice acceptor
        self.variant["Pos"] = "431254483"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceAcceptor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_getVarSpliceRegionBoundsAcceptorBRCA1(self, varInSpliceRegion, getRefSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position) 
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        donor = False

        # checks that variant in exon 21 splice acceptor region boundaries are returned correctly
        self.variant["Pos"] = "43051117"
        spliceAccRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor)
        self.assertEquals(exonAcceptorBoundsBRCA1["exon21"]["acceptorStart"], spliceAccRegion["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA1["exon21"]["acceptorEnd"], spliceAccRegion["acceptorEnd"])

    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_varInSpliceRegionAcceptorBRCA2(self, getRefSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Variant is correctly identified as in or NOT in a splice acceptor region
           uses multiple positions across multiple exons for BRCA2
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = False

        # checks that -21st base in intron is NOT counted as in splice acceptor
        self.variant["Pos"] = "32357721"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceAcceptor)

        # checks that first base in intron counted as in splice acceptor
        self.variant["Pos"] = "32316402"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in intron counted as in splice acceptor
        self.variant["Pos"] = "32319069"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in intron counted as in splice acceptor
        self.variant["Pos"] = "32325075"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that first base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326101"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326243"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326501"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertTrue(inSpliceAcceptor)

        # checks that 4th base in exon is NOT counted as in splice acceptor
        self.variant["Pos"] = "32362526"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceAcceptor)

        # checks that region before exon 1 is NOT counted as in splice acceptor
        self.variant["Pos"] = "32315479"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor)
        self.assertFalse(inSpliceAcceptor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_getVarSpliceRegionBoundsAcceptorBRCA2(self, varInSpliceRegion, getRefSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position) 
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = False

        # checks that variant in exon 20 splice acceptor region boundaries are returned correctly
        self.variant["Pos"] = "32370948"
        spliceAccRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor)
        self.assertEquals(exonAcceptorBoundsBRCA2["exon20"]["acceptorStart"], spliceAccRegion["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA2["exon20"]["acceptorEnd"], spliceAccRegion["acceptorEnd"])

    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCiDomainEnigmaBRCA1(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA1 as defined by ENIGMA rules'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks variant in BRCA 1 RING domain is identified as in ENIGMA CI domain
        self.variant["Pos"] = "43124089"
        inEnigmaCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inEnigmaCI)

        # checks variant in BRCA1 BRCT domain is identified as in ENIGMA CI domain
        self.variant["Pos"] = "43070945"
        inEnigmaCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inEnigmaCI)

        # checks variant NOT in BRCA1 CI domain is NOT identified as in ENIGMA CI domain
        self.variant["Pos"] = "43097274"
        inEnigmaCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertFalse(inEnigmaCI)
        
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCiDomainEnigmaBRCA2(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA2 as defined by ENIGMA rules'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks variant in BRCA2 DNB domain is identified as in ENIGMA CI domain
        self.variant["Pos"] = "32379809"
        inEnigmaCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inEnigmaCI)

        # checks variant NOT in BRCA2 CI domain is NOT identified as in ENIGMA CI domain
        self.variant["Pos"] = "32398354"
        inEnigmaCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertFalse(inEnigmaCI)

    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCiDomainPriorsBRCA1(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA1 as defined by PRIORS webiste'''
        boundaries = "priors"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks variant in BRCA1 initiation domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "43124096"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA1 RING domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "43115746"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA1 BRCT domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "43057092"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks that variant NOT in BRCA1 CI domain is NOT identified as in PRIORS CI domain
        self.variant["Pos"] = "43124090"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertFalse(inPriorsCI)
        
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCiDomainPriorsBRCA2(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA2 as defined by PRIORS webiste'''
        boundaries = "priors"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks variant in BRCA2 initiation domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32316462"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA2 PALB2 domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32319092"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA2 DNB domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32362561"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in TR2/RAD5 domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32398406"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks that variant NOT in BRCA2 CI domain is NOT identified as in PRIORS CI domain
        self.variant["Pos"] = "32336283"
        inPriorsCI = calcVarPriors.varInCiDomain(self.variant, boundaries)
        self.assertFalse(inPriorsCI)

    def test_varInGreyZone(self):
        '''Tests that variant is correctly identified as in the grey zone'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks that BRCA1 variant is NOT considered in grey zone
        self.variant["Pos"] = "43045708"
        inGreyZone = calcVarPriors.varInGreyZone(self.variant)
        self.assertFalse(inGreyZone)

        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks that variant before the BRCA2 grey zone is NOT identified as in BRCA2 grey zone
        self.variant["Pos"] = "32398437"
        inGreyZone = calcVarPriors.varInGreyZone(self.variant)
        self.assertFalse(inGreyZone)

        # checks that variant in BRCA2 grey zone is identified as in BRCA2 grey zone
        self.variant["Pos"] = "32398459"
        inGreyZone = calcVarPriors.varInGreyZone(self.variant)
        self.assertTrue(inGreyZone)
        
    def test_varAfterGreyZoneBRCA1(self):
        '''Tests that variant in BRCA1 is NOT considered as after the grey zone'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks that BRCA1 variant is NOT considered after grey zone
        self.variant["Pos"] = "43045689"
        afterGreyZone = calcVarPriors.varAfterGreyZone(self.variant)
        self.assertFalse(afterGreyZone)

    @mock.patch('calcVarPriors.varInUTR', return_value = False)
    @mock.patch('calcVarPriors.varInGreyZone', return_value = True)
    def test_varAfterGreyZoneFalseBRCA2(self, varInUTR, varInGreyZone):
        '''Tests that variant in BRCA2 is correctly identified as NOT after the grey zone'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        self.variant["Pos"] = "32398459"
        afterGreyZone = calcVarPriors.varAfterGreyZone(self.variant)
        self.assertFalse(afterGreyZone)
        
    @mock.patch('calcVarPriors.varInUTR', return_value = False)
    @mock.patch('calcVarPriors.varInGreyZone', return_value = False)
    def test_varAfterGreyZoneFalseBRCA2(self, varInUTR, varInGreyZone):
        '''Tests that variant after BRCA2 grey zone is correctly identified as after the grey zone'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        
        self.variant["Pos"] = "32398489"
        afterGreyZone = calcVarPriors.varAfterGreyZone(self.variant)
        self.assertTrue(afterGreyZone)

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = True)
    def test_getVarLocationOutBounds(self, varOutsideBoundaries):
        '''Tests that variants outside transcript boundaries are correctly identified as outside transcript boundaries'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        boundaries = "enigma"
        
        # BRCA1 variant outside transcript boundaries (before txn start)
        self.variant["Pos"] = "43125600"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["outBounds"])

        # BRCA1 variant outside transcript boundaries (after txn end)
        self.variant["Pos"] = "43044000"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["outBounds"])

        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"
        
        # BRCA2 variant outside transcript boundaries (before txn start)
        self.variant["Pos"] = "32315300"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["outBounds"])

        # BRCA2 variant outside transcript boundaries (after txn end)
        self.variant["Pos"] = "32399800"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["outBounds"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInCiDomain', return_value = True)
    def test_getVarLocationCiDomain(self, varOutsideBoundaries, varInExon, varInSpliceRegion, varInCiDomain):
        '''
        Tests that variants in either PRIORS or ENIGMA CI domain and NOT in a splice region are
        identified as in CI domain
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        boundaries = "enigma"
        
        # BRCA1 variant in middle of ENIGMA CI domain
        self.variant["Pos"] = "43063930"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCI"])

        # BRCA1 variant in middle of PRIORS CI domain
        boundaries = "priors"
        self.variant["Pos"] = "43106502"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCI"])
        
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"
        
        # BRCA2 variant in middle of ENIGMA CI domain
        self.variant["Pos"] = "32376714"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCI"])

        # BRCA2 variant in middle of PRIORS CI domain
        boundaries = "priors"
        self.variant["Pos"] = "32363207"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCI"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInCiDomain', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca1RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_getVarLocationCiDomainSpliceRegionBRCA1(self, varOutsideBoundaries, varInExon, varInCiDomain,
                                                     getRefSpliceDonorBoundaries, getRefSpliceAcceptorBoundaries):
        '''
        Tests that BRCA1 variants in either PRIORS or ENIGMA CI domains AND in splice region are
        correctly identified as in CI domain splice region
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
                
        # BRCA1 variant in ENIGMA CI domain splice donor region
        boundaries = "enigma"
        self.variant["Pos"] = "43104868"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceDonor"])
        
        # BRCA1 variant in PRIORS CI domain splice donor region
        boundaries = "priors"
        self.variant["Pos"] = "43106456"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceDonor"])

        # BRCA1 variant in ENIGMA CI domain splice acceptor region
        boundaries = "enigma"
        self.variant["Pos"] = "43063373"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceAcceptor"])

        # BRCA1 variant in PRIORS CI domain splice acceptor region
        boundaries = "priors"
        self.variant["Pos"] = "43057135"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceAcceptor"])
        
    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInCiDomain', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_getVarLocationCiDomainSpliceRegionBRCA2(self, varOutsideBoundaries, varInExon, varInCiDomain,
                                                     getRefSpliceDonorBoundaries, getRefSpliceAcceptorBoundaries):
        '''
        Tests that BRCA2 variants in either PRIORS or ENIGMA CI domains AND in splice region are
        correctly identified as in CI domain splice region
        '''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # BRCA2 variant in ENIGMA CI splice donor region
        boundaries = "enigma"
        self.variant["Pos"] = "32357929"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceDonor"])

        # BRCA2 variant in PRIORS CI splice donor region
        boundaries = "priors"
        self.variant["Pos"] = "32316527"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceDonor"])

        # BRCA2 variant in ENIGMA CI splice acceptor region
        boundaries = "enigma"
        self.variant["Pos"] = "32376670"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceAcceptor"])

        # BRCA2 variant in PRIORS CI splice acceptor region
        boundaries = "priors"
        self.variant["Pos"] = "32363179"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inCISpliceAcceptor"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca1RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_getVarLocationSpliceRegionBRCA1(self, varOutsideBoundaries, getExonBoundaries, getRefSpliceDonorBoundaries,
                                             getRefSpliceAcceptorBoundaries):
        '''Tests that BRCA1 variants in splice regions are correctly identified as in splice regions'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        boundaries = "enigma"
        
        # BRCA1 variant in splice donor region
        self.variant["Pos"] = "43074331"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inSpliceDonor"])
        
        # BRCA1 variant in splice acceptor region
        self.variant["Pos"] = "43082575"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inSpliceAcceptor"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.getRefSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_getVarLocationSpliceRegionBRCA2(self, varOutsideBoundaries, getExonBoundaries, getRefSpliceDonorBoundaries,
                                             getRefSpliceAcceptorBoundaries):
        '''Tests that BRCA1 variants in splice regions are correctly identified as in splice regions'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"
        
        # BRCA2 variant in splice donor region
        self.variant["Pos"] = "32333388"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inSpliceDonor"])
        
        # BRCA2 variant in splice acceptor region
        self.variant["Pos"] = "32329443"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inSpliceAcceptor"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInCiDomain', return_value = False)
    def test_getVarLocationInExon(self, varOutsideBoundaries, varInExon, varInSpliceRegion, varInCiDomain):
        '''Tests that variants in exons are correctly identified as in exons'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        boundaries = "enigma"

        # BRCA1 variant in exon
        self.variant["Pos"] = "43071220"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inExon"])
        
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"

        # BRCA2 variant in exon
        self.variant["Pos"] = "32336289"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inExon"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInUTR', return_value = True)
    def test_getVarLocationInUtrBRCA1(self, varOutsideBoundaries, getExonBoundaries, varInSpliceRegion, varInUTR):
        '''Tests that BRCA1 variants in UTR are correctly identified as in UTR'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        boundaries = "enigma"

        # BRCA1 variant in 5' UTR
        self.variant["Pos"] = "43124138"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inUTR"])

        # BRCA1 variant in 3' UTR
        self.variant["Pos"] = "43045660"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inUTR"])
        
    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInUTR', return_value = True)
    def test_getVarLocationInUtrBRCA2(self, varOutsideBoundaries, getExonBoundaries, varInSpliceRegion, varInUTR):
        '''Tests that BRCA2 variants in UTR are correctly identified as in UTR'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"

        # BRCA2 variant in 5' UTR
        self.variant["Pos"] = "32316398"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inUTR"])

        # BRCA2 variant in 3' UTR
        self.variant["Pos"] = "32398790"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inUTR"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    def test_getVarLocationInIntron(self, varOutsideBoundaries, varInExon, varInSpliceRegion):
        '''Tests that variants in introns are correctly identified as in introns'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        boundaries = "enigma"

        # BRCA1 variant in intron
        self.variant["Pos"] = "43071263"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inIntron"])

        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"
        
        # BRCA2 variant in intron
        self.variant["Pos"] = "32344537"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inIntron"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInGreyZone', return_value = True)
    def test_getVarLocationInGreyZone(self, varOutsideBoundaries, varInExon, varInSpliceRegion, varInGreyZone):
        '''Tests that BRCA2 variant in grey zone is correctly identified as in grey zone'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"

        # BRCA2 variant in grey zone
        self.variant["Pos"] = "32398465"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["inGreyZone"])

    @mock.patch('calcVarPriors.varOutsideBoundaries', return_value = False)
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varAfterGreyZone', return_value = True)
    def test_getVarLocationAfterGreyZone(self, varOutsideBoundaries, varInExon, varInSpliceRegion, varAfterGreyZone):
        '''Tests that BRCA2 variant after grey zone is correctly identified as after grey zone'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        boundaries = "enigma"

        # BRCA2 variant after grey zone
        self.variant["Pos"] = "32398499"
        varLoc = calcVarPriors.getVarLocation(self.variant, boundaries)
        self.assertEquals(varLoc, variantLocations["afterGreyZone"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca1Seq)    
    def test_getSeqLocDictBRCA1(self, getFastaSeq):
        '''
        Tests that boundary endpoints for BRCA1 are:
        1. Included in dictionary
        2. Have the correct base based on mocked sequence
        '''
        chrom = "chr17"
        strand = "-"
        rangeStart = 43051137
        rangeStop = 43051115
        seqLocDict = calcVarPriors.getSeqLocDict(chrom, strand, rangeStart, rangeStop)
        self.assertEquals(seqLocDict[rangeStart], brca1Seq[-1])
        self.assertEquals(seqLocDict[rangeStop], brca1Seq[0])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca2Seq)
    def test_getSeqLocDictBRCA2(self, getFastaSeq):
        '''
        Tests that boundary endpoints for BRCA2 are:
        1. Included in dictionary
        2. Have the correct base based on mocked sequence
        '''
        chrom = "chr13"
        strand = "+"
        rangeStart = 32370936
        rangeStop = 32370958
        seqLocDict = calcVarPriors.getSeqLocDict(chrom, strand, rangeStart, rangeStop)
        self.assertEquals(seqLocDict[rangeStart], brca2Seq[0])
        self.assertEquals(seqLocDict[rangeStop], brca2Seq[-1])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca1Seq)
    def test_getAltSeqDictBRCA1(self, getFastaSeq):
        '''
        Tests that for given variant genomic position:
        1. Reference allele is set correctly in original dictionary (refSeqDict)
        2. Alternate allele is set correctly in alternate sequence dicitonary (altSeqDict)
        '''
        chrom = "chr17"
        strand = "-"
        rangeStart = 43051137
        rangeStop = 43051115
        refSeqDict = calcVarPriors.getSeqLocDict(chrom, strand, rangeStart, rangeStop)
        self.variant["Pos"] = "43051120"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        altSeqDict = calcVarPriors.getAltSeqDict(self.variant, refSeqDict)
        self.assertEquals(refSeqDict[int(self.variant["Pos"])], self.variant["Ref"])
        self.assertEquals(altSeqDict[int(self.variant["Pos"])], self.variant["Alt"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca2Seq)
    def test_getAltSeqDictBRCA2(self, getFastaSeq):
        '''
        Tests that for given variant genomic position:
        1. Reference allele is set correctly in original dictionary (refSeqDict)
        2. Alternate allele is set correctly in alternate sequence dictionary (altSeqDict)
        '''
        chrom = "chr13"
        strand = "+"
        rangeStart = 32370936
        rangeStop = 32370958
        refSeqDict = calcVarPriors.getSeqLocDict(chrom, strand, rangeStart, rangeStop)
        self.variant["Pos"] = "32370944"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        altSeqDict = calcVarPriors.getAltSeqDict(self.variant, refSeqDict)
        self.assertEquals(refSeqDict[int(self.variant["Pos"])], self.variant["Ref"])
        self.assertEquals(altSeqDict[int(self.variant["Pos"])], self.variant["Alt"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca1Seq)
    def test_getAltSeqBRCA1(self, getFastaSeq):
        '''Tests that alternate sequence string generated is correct for - strand gene (BRCA1)'''
        chrom = "chr17"
        strand = "-"
        rangeStart = 43051137
        rangeStop = 43051115
        refSeqDict = calcVarPriors.getSeqLocDict(chrom, strand, rangeStart, rangeStop)
        self.variant["Pos"] = "43051120"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        altSeqDict = calcVarPriors.getAltSeqDict(self.variant, refSeqDict)
        altSeq = calcVarPriors.getAltSeq(altSeqDict, strand)
        # reference sequence on plus strand
        brca1RefSeq = "GATCTGGAAGAAGAGAGGAAGAG"
        # reverse complement with variant included
        brca1AltSeq = "CTCTTCCTCTCTTCTTCGAGATC"
        self.assertEquals(altSeq, brca1AltSeq)

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca2Seq)
    def test_getAltSeqBRCA2(self, getFastaSeq):
        '''Tests that alternate sequence string generated is correct for + strand gene (BRCA2)'''
        chrom = "chr13"
        strand = "+"
        rangeStart = 32370936
        rangeStop = 32370958
        refSeqDict = calcVarPriors.getSeqLocDict(chrom, strand, rangeStart, rangeStop)
        self.variant["Pos"] = "32370944"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        altSeqDict = calcVarPriors.getAltSeqDict(self.variant, refSeqDict)
        altSeq = calcVarPriors.getAltSeq(altSeqDict, strand)
        # reference sequence on plus strand
        brca2RefSeq = "TGTGTAACACATTATTACAGTGG"
        # alternate sequence containng the alterante allele
        brca2AltSeq = "TGTGTAACCCATTATTACAGTGG"
        self.assertEquals(altSeq, brca2AltSeq)

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca1Seq)
    def test_getRefAltSeqsBRCA1(self, getFastaSeq):
        '''Tests that ref and alt sequence are generated correctly for - strand gene (BRCA1)'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        rangeStart = 43051137
        rangeStop = 43051115
        self.variant["Pos"] = "43051120"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        refAltSeqs = calcVarPriors.getRefAltSeqs(self.variant, rangeStart, rangeStop)
        # reference sequence on minus strand
        brca1RefSeq = "CTCTTCCTCTCTTCTTCCAGATC"
        # reverse complement with variant included
        brca1AltSeq = "CTCTTCCTCTCTTCTTCGAGATC"
        self.assertEquals(refAltSeqs["refSeq"], brca1RefSeq)
        self.assertEquals(refAltSeqs["altSeq"], brca1AltSeq)

    @mock.patch('calcVarPriors.getFastaSeq', return_value = brca2Seq)
    def test_getRefAltSeqsBRCA2(self, getFastaSeq):
        '''Tests that ref and alt sequence are generated correctly for + strand gene (BRCA2)'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        rangeStart = 32370936
        rangeStop = 32370958
        self.variant["Pos"] = "32370944"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        refAltSeqs = calcVarPriors.getRefAltSeqs(self.variant, rangeStart, rangeStop)
        # reference sequence on plus strand
        brca2RefSeq = "TGTGTAACACATTATTACAGTGG"
        # alternate sequence containng the alterante allele
        brca2AltSeq = "TGTGTAACCCATTATTACAGTGG"
        self.assertEquals(refAltSeqs["refSeq"], brca2RefSeq)
        self.assertEquals(refAltSeqs["altSeq"], brca2AltSeq)

    def test_getZScore(self):
        '''
        Tests that:
        1. For score in splice donor site:
            - checks that zscore is less than zero if MaxEntScan score less than donor mean
            - checks that zscore is greater than zero if MaxEntScan score is greater than donor mean
        2. For score in splice acceptor site:
            - checks that zscore is less than zero if MaxEntScan score is less than acceptor mean
            - checks that zscore is greater than zero if MaxEntScan score is greater than acceptor mean
        '''
        # score less than donor mean of ~7.94
        maxEntScanScore = 7.8
        zScore = calcVarPriors.getZScore(maxEntScanScore, donor=True)
        self.assertLess(zScore, 0)

        # score greater than donor mean of ~7.94
        maxEntScanScore = 7.99
        zScore = calcVarPriors.getZScore(maxEntScanScore, donor=True)
        self.assertGreater(zScore, 0)

        # score less than acceptor mean of ~7.98
        maxEntScanScore = 7.9
        zScore = calcVarPriors.getZScore(maxEntScanScore, donor=False)
        self.assertLess(zScore, 0)

        # score less than acceptor mean of ~7.98
        maxEntScanScore = 8
        zScore = calcVarPriors.getZScore(maxEntScanScore, donor=False)
        self.assertGreater(zScore, 0)
        
    def test_getEnigmaClass(self):
        ''''
        Tests that predicted qualititative ENIGMA class is assigned correctly based on prior prob
        Specifically tests for priors in class 1 and class 5
        and most commonly assigned priorProb = 0.04, 0.34, and 0.97
        '''
        priorProb = 0.0001
        enigmaClass = calcVarPriors.getEnigmaClass(priorProb)
        self.assertEquals(enigmaClass, enigmaClasses["class1"])
        
        priorProb = 0.04
        enigmaClass = calcVarPriors.getEnigmaClass(priorProb)
        self.assertEquals(enigmaClass, enigmaClasses["class2"])

        priorProb = 0.34
        enigmaClass = calcVarPriors.getEnigmaClass(priorProb)
        self.assertEquals(enigmaClass, enigmaClasses["class3"])

        priorProb = 0.97
        enigmaClass = calcVarPriors.getEnigmaClass(priorProb)
        self.assertEquals(enigmaClass, enigmaClasses["class4"])
                
        priorProb = 0.995
        enigmaClass = calcVarPriors.getEnigmaClass(priorProb)
        self.assertEquals(enigmaClass, enigmaClasses["class5"])

        priorProb = "N/A"
        enigmaClass = calcVarPriors.getEnigmaClass(priorProb)
        self.assertEquals(enigmaClass, None)

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCTTACCTT")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceDonorBounds["exon14"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.57,
                                                                               "zScore": 1.1300618149879533},
                                                                 "altScores": {"maxEntScanScore": 9.21,
                                                                               "zScore": 0.5461191272666394}})
    def test_getPriorProbSpliceDonorSNSLowProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                    getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variant that creates a resaonble splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that creates a reasonable splice donor site
        self.variant["Pos"] = "43076489"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCTTACCTT")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceDonorBounds["exon14"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.57,
                                                                               "zScore": 1.1300618149879533},
                                                                 "altScores": {"maxEntScanScore": 6.12,
                                                                               "zScore": -0.7806330088060529}})
    def test_getPriorProbSpliceDonorSNSModerateProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                         getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variant that weakens a reasonably strong splice donor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that weakens a reasonably strong splice donor site
        self.variant["Pos"] = "43076485"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["moderate"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TTTTACCAA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceDonorBounds["exon7"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 3.23,
                                                                               "zScore": -2.021511220213846},
                                                                 "altScores": {"maxEntScanScore": -4.42,
                                                                               "zScore": -5.306188838646238}})
    def test_getPriorProbSpliceDonorSNSHighProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                     getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA1 variant that further weakens a weak splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that further weakens a weak splice donor site
        self.variant["Pos"] = "43104120"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["high"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class4"])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCTTACCTT")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceDonorBounds["exon14"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.57,
                                                                               "zScore": 1.1300618149879533},
                                                                 "altScores": {"maxEntScanScore": 10.77,
                                                                               "zScore": 1.2159357396528523}})
    def test_getPriorProbSpliceDonorSNSImprovedProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                         getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variant that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that makes a splice donor site stronger or equally strong
        self.variant["Pos"] = "43076490"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCGGTAAGA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceDonorBounds["exon13"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.53,
                                                                               "zScore": 1.1128870300549731},
                                                                 "altScores": {"maxEntScanScore": 8.91,
                                                                               "zScore": 0.4173082402692903}})
    def test_getPriorProbSpliceDonorSNSLowProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                    getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that creates a resaonble splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that creates a reasonable splice donor site
        self.variant["Pos"] = "32346895"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCGGTAAGA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceDonorBounds["exon13"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.53,
                                                                               "zScore": 1.1128870300549731},
                                                                 "altScores": {"maxEntScanScore": 4.35,
                                                                               "zScore": -1.5406172420904107}})
    def test_getPriorProbSpliceDonorSNSModerateProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                         getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that weakens a reasonably strong splice donor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that weakens a reasonably strong splice donor site
        self.variant["Pos"] = "32346899"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["moderate"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "CAGGCAAGT")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceDonorBounds["exon17"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 3.1,
                                                                               "zScore": -2.07732927124603},
                                                                 "altScores": {"maxEntScanScore": 0.56,
                                                                               "zScore": -3.1679281144902496}})
    def test_getPriorProbSpliceDonorSNSHighProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                     getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA2 variant that further weakens a weak splice donor site'''  
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that further weakens a weak splice donor site
        self.variant["Pos"] = "32362693"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["high"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class4"])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCGGTAAGA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_donor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceDonorBounds["exon13"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.53,
                                                                               "zScore": 1.1128870300549731},
                                                                 "altScores": {"maxEntScanScore": 11.78,
                                                                               "zScore": 1.649599059210593}})
    def test_getPriorProbSpliceDonorSNSImprovedProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                         getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that makes a splice donor site stronger or equally strong
        self.variant["Pos"] = "32346902"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "CATCTGTAAAATACAAGGGAAAA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceAcceptorBounds["exon7"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 11.68,
                                                                               "zScore": 1.5183252360035546},
                                                                 "altScores": {"maxEntScanScore": 10.94,
                                                                               "zScore": 1.214256756422072}})
    def test_getPriorProbSpliceAcceptorSNSLowProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                       getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variants that creates a resaonble splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that creates a reasonable splice acceptor site
        self.variant["Pos"] = "43104275"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "CATCTGTAAAATACAAGGGAAAA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceAcceptorBounds["exon7"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 11.68,
                                                                               "zScore": 1.5183252360035546},
                                                                 "altScores": {"maxEntScanScore": 9.01,
                                                                               "zScore": 0.4212132894055031}})
    def test_getPriorProbSpliceAcceptorSNSModerateProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variants that weakens a reasonably strong splice acceptor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that weakens a reasonably strong splice acceptor site
        self.variant["Pos"] = "43104267"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["moderate"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "GAACTTTAACACATTAGAAAAAC")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceAcceptorBounds["exon2"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 4.9,
                                                                               "zScore":  -1.2675994823240817},
                                                                 "altScores": {"maxEntScanScore": -3.17,
                                                                               "zScore": -4.583589523165384}})
    def test_getPriorProbSpliceAcceptorSNSHighProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                        getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA1 variant that further weakens a weak splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that further weakens a weak splice acceptor site
        self.variant["Pos"] = "43124116"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["high"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class4"])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "CATCTGTAAAATACAAGGGAAAA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca1RefSpliceAcceptorBounds["exon7"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 11.68,
                                                                               "zScore": 1.5183252360035546},
                                                                 "altScores": {"maxEntScanScore": 12.41,
                                                                               "zScore": 1.8182846820771794}})
    def test_getPriorProbSpliceAcceptorSNSImprovedProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variants that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
    
        # checks prior prob for BRCA1 variant that makes a splice acceptor site stronger or equally strong
        self.variant["Pos"] = "43104261"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCTCATCTTTCTCCAAACAGTTA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceAcceptorBounds["exon23"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.35,
                                                                               "zScore": 0.9718237794584578},
                                                                 "altScores": {"maxEntScanScore": 10.09,
                                                                               "zScore": 0.8649889082541532}})
    def test_getPriorProbSpliceAcceptorSNSLowProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                       getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that creates a resaonble splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that creates a reasonable splice acceptor site
        self.variant["Pos"] = "32379751"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCTCATCTTTCTCCAAACAGTTA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceAcceptorBounds["exon23"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.35,
                                                                               "zScore": 0.9718237794584578},
                                                                 "altScores": {"maxEntScanScore": 8.84,
                                                                               "zScore": 0.3513597197719193}})
    def test_getPriorProbSpliceAcceptorSNSModerateProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that weakens a reasonably strong splice acceptor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that weakens a reasonably strong splice acceptor site
        self.variant["Pos"] = "32379746"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["moderate"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getFastaSeq', return_value = "AAGTATTTATTCTTTGATAGATT")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceAcceptorBounds["exon15"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 5.16,
                                                                               "zScore": -1.1607646111197771},
                                                                 "altScores": {"maxEntScanScore": -2.91,
                                                                               "zScore": -4.4767546519610795}})
    def test_getPriorProbSpliceAcceptorSNSHighProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                        getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA2 variant that further weakens a weak splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that further weakens a weak splice acceptor site
        self.variant["Pos"] = "32356427"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["high"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class4"])
        
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "TCTCATCTTTCTCCAAACAGTTA")
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = "splice_acceptor_variant")
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = brca2RefSpliceAcceptorBounds["exon23"])
    @mock.patch('calcVarPriors.getRefAltScores', return_value = {"refScores": {"maxEntScanScore": 10.35,
                                                                               "zScore": 0.9718237794584578},
                                                                 "altScores": {"maxEntScanScore": 10.42,
                                                                               "zScore": 1.000587014013463}})
    def test_getPriorProbSpliceAcceptorSNSImprovedProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that makes a splice acceptor site stronger or equally strong
        self.variant["Pos"] = "32379731"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["low"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = variantLocations["afterGreyZone"])
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "missense_variant")
    def test_getPriorProbAfterGreyZoneMissenseSNS(self, getVarType, getVarLocation, getVarConsequences):
        '''
        Tests that:
        prior prob is set to N/A and ENIGMA class is class 2 for a BRCA2 missense variant after the grey zone
        '''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Chr"] = "13"
        self.variant["Hg38_Start"] = "32398528"
        self.variant["Hg38_End"] = "32398528"

        self.variant["Pos"] = "32398528"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbAfterGreyZoneSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.getVarLocation', return_value = variantLocations["afterGreyZone"])
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    def test_getPriorProbAfterGreyZoneNonesenseSNS(self, getVarType, getVarLocation, getVarConsequences):
        '''
        Tests that:
        prior prob is set to N/A and ENIGMA class is class 2 for a BRCA2 nonsense variant after the grey zone
        '''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Chr"] = "13"
        self.variant["Hg38_Start"] = "32398492"
        self.variant["Hg38_End"] = "32398492"

        self.variant["Pos"] = "32398492"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbAfterGreyZoneSNS(self.variant, boundaries)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcMaxEntScanMeanStd.fetch_gene_coordinates', return_value = transcriptDataBRCA2)
    def test_getVarDict(self, fetch_gene_coordinates):
        '''
        Tests that: 
        1. Variant information is being parsed correctly
        '''
        boundaries = "enigma"
        varDict = calcVarPriors.getVarDict(self.variant, boundaries)
        self.assertEquals(varDict["varHGVScDNA"], self.variant["pyhgvs_cDNA"])
        self.assertEquals(varDict["varChrom"], self.variant["Chr"])
        self.assertEquals(varDict["varGene"], self.variant["Gene_Symbol"])
        self.assertEquals(varDict["varGenCoordinate"], self.variant["Pos"])
        

