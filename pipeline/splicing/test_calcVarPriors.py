import unittest
import mock
import calcVarPriors
from calcVarPriors import STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH
from calcVarPriors import STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH, STD_DE_NOVO_OFFSET
import calcMaxEntScanMeanStd
from calcVarPriorsMockedResponses import brca1Exons, brca2Exons 
from calcVarPriorsMockedResponses import brca1RefSpliceDonorBounds, brca2RefSpliceDonorBounds 
from calcVarPriorsMockedResponses import brca1RefSpliceAcceptorBounds, brca2RefSpliceAcceptorBounds
from calcVarPriorsMockedResponses import variantData

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
                 "class5": "class_5",
                 "NA": "N/A"}

# possible prior probability of pathogenecity values
priorProbs = {"deNovoLow": 0.02,
              "proteinLow": 0.03,
              "low": 0.04,
              "proteinMod": 0.29,
              "deNovoMod": 0.3,
              "moderate": 0.34,
              "proteinHigh": 0.81,
              "deNovoHigh": 0.64,
              "high": 0.97,
              "pathogenic": 0.99,
              "NA": "N/A"}


class test_calcVarPriors(unittest.TestCase):

    def setUp(self):

        self.variant = {"Chr": "13",
                        "Pos": "32314943",
                        "Ref": "A",
                        "Alt": "G",
                        "Gene_Symbol": "BRCA2",
                        "Reference_Sequence": "NM_000059.3",
                        "HGVS_cDNA": "c.-764A>G"}
                  
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
        Uses example boundaries defined at beginning of script
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        spliceDonorBounds = calcVarPriors.getRefSpliceDonorBoundaries(self.variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
        # checks that region after last exon is not considered a splice donor region
        self.assertNotIn("exon24", spliceDonorBounds)
        # to find exon specified in global variables
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
        Uses example boundaries defined at beginning of script
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        spliceDonorBounds = calcVarPriors.getRefSpliceDonorBoundaries(self.variant, STD_DONOR_INTRONIC_LENGTH, STD_DONOR_EXONIC_LENGTH)
        # checks that region after last exon is not considered a splice donor region
        self.assertNotIn("exon27", spliceDonorBounds)
        # to find exon specified in global variables
        exon = exonDonorBoundsBRCA2.keys()[0]
        self.assertEquals(exonDonorBoundsBRCA2[exon]["donorStart"],
                          spliceDonorBounds[exon]["donorStart"])
        self.assertEquals(exonDonorBoundsBRCA2[exon]["donorEnd"],
                          spliceDonorBounds[exon]["donorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_getSpliceAcceptorBoundariesRefBRCA1(self, getExonBoundaries, getVarStrand):
        '''
        Tests that ref splice acceptor boundaries are set correctly for reference transcript (NM_007294.3) and strand (-)
        Uses example boundaries defined at beginning of script
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        spliceAcceptorBounds = calcVarPriors.getSpliceAcceptorBoundaries(self.variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
        # checks that region before first exon is not considered a splice acceptor region
        self.assertNotIn("exon1", spliceAcceptorBounds)
        # to find exon specified in global variables
        exon = exonAcceptorBoundsBRCA1.keys()[0]
        self.assertEquals(exonAcceptorBoundsBRCA1[exon]["acceptorStart"],
                          spliceAcceptorBounds[exon]["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA1[exon]["acceptorEnd"],
                          spliceAcceptorBounds[exon]["acceptorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_getSpliceAcceptorBoundariesDeNovoBRCA1(self, getExonBoundaries, getVarStrand):
        '''
        Tests that de novo splice acceptor boundaries are set correctly for reference transcript (NM_007294.3) and strand (-)
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        deNovoSpliceAccBounds = calcVarPriors.getSpliceAcceptorBoundaries(self.variant, STD_ACC_INTRONIC_LENGTH, STD_DE_NOVO_LENGTH)
        expectedDeNovoRegionExon6 = {"acceptorStart": 43104976,
                                     "acceptorEnd": 43104947}
        self.assertEquals(deNovoSpliceAccBounds["exon6"]["acceptorStart"],
                          expectedDeNovoRegionExon6["acceptorStart"])
        self.assertEquals(deNovoSpliceAccBounds["exon6"]["acceptorEnd"],
                          expectedDeNovoRegionExon6["acceptorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_getSpliceAcceptorBoundariesRefBRCA2(self, getExonBoundaries, getVarStrand):
        '''
        Tests that ref splice acceptor boundaries are set correctly for reference transcript (NM_000059.3) and strand (+)
        Uses example boundaries defined at beginning of script
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        spliceAcceptorBounds = calcVarPriors.getSpliceAcceptorBoundaries(self.variant, STD_ACC_INTRONIC_LENGTH, STD_ACC_EXONIC_LENGTH)
        # checks that region before first exon is not considered a splice acceptor region
        self.assertNotIn("exon1", spliceAcceptorBounds)
        # to find exon specified in global variables
        exon = exonAcceptorBoundsBRCA2.keys()[0]
        self.assertEquals(exonAcceptorBoundsBRCA2[exon]["acceptorStart"],
                          spliceAcceptorBounds[exon]["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA2[exon]["acceptorEnd"],
                          spliceAcceptorBounds[exon]["acceptorEnd"])

    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_getSpliceAcceptorBoundariesDeNovoBRCA2(self, getExonBoundaries, getVarStrand):
        '''
        Tests that de novo splice acceptor boundaries are set correctly for reference transcript (NM_000059.3) and strand (+)
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        deNovoSpliceAccBounds = calcVarPriors.getSpliceAcceptorBoundaries(self.variant, STD_ACC_INTRONIC_LENGTH, STD_DE_NOVO_LENGTH)
        expectedDeNovoRegionExon8 = {"acceptorStart": 32329423,
                                     "acceptorEnd": 32329452}
        self.assertEquals(deNovoSpliceAccBounds["exon8"]["acceptorStart"],
                          expectedDeNovoRegionExon8["acceptorStart"])
        self.assertEquals(deNovoSpliceAccBounds["exon8"]["acceptorEnd"],
                          expectedDeNovoRegionExon8["acceptorEnd"])
        
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

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_getVarExonNumberSNSBRCA1(self, varInExon, getExonBoundaries, getVarStrand):
        '''Tests that exon number is set correctly for minus strand (BRCA1) variant in exon'''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"

        # variant position in exon 13
        self.variant["Pos"] = "43082564"
        print self.variant
        varExonNum = calcVarPriors.getVarExonNumberSNS(self.variant)
        self.assertEquals(varExonNum, "exon13")

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_getVarExonNumberSNSBRCA2(self, varInExon, getExonBoundaries, getVarStrand):
        '''Tests that exon number is set correctly for plus strand (BRCA2) variant in exon'''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"

        # variant position in exon 4
        self.variant["Pos"] = "32325166"
        varExonNum = calcVarPriors.getVarExonNumberSNS(self.variant)
        self.assertEquals(varExonNum, "exon4")
        
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
        deNovo = False

        #checks that 7th base in intron is NOT counted as in splice donor
        self.variant["Pos"] = "43063326"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceDonor)
        
        # checks that first base in intron counted as in splice donor
        self.variant["Pos"] =  "43097243"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in intron counted as in splice donor
        self.variant["Pos"] = "43095842"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that last base in intron counted as in splice donor
        self.variant["Pos"] = "43091429"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that first base in exon counted as in splice donor
        self.variant["Pos"] = "43090946"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in exon counted as in splice donor
        self.variant["Pos"] = "43082405"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that last base in exon counted as in splice donor
        self.variant["Pos"] = "43076488"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that 4th to last base in exon is NOT counted as in splice donor
        self.variant["Pos"] = "43057055"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceDonor)

        # checks that region after exon 24 is  NOT counted as in splice donor
        self.variant["Pos"] = "43044294"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
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
        deNovo = False
    
        # checks that variant in exon 16 splice donor region boundaries are returned correctly
        self.variant["Pos"] = "43070925"
        spliceDonorRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor, deNovo=deNovo)
        self.assertEquals(spliceDonorRegion["exonName"], "exon16")
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
        deNovo = False

        # checks that 7th base in intron is NOT counted as in splice donor
        self.variant["Pos"] = "32363540"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceDonor)

        # checks that first base in intron counted as in splice donor
        self.variant["Pos"] = "32329493"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in intron counted as in splice donor
        self.variant["Pos"] = "32331033"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that last base in intron counted as in splice donor
        self.variant["Pos"] = "32333393"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that first base in exon counted as in splice donor
        self.variant["Pos"] = "32341194"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in exon counted as in splice donor
        self.variant["Pos"] = "32344652"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that last base in exon counted as in splice donor
        self.variant["Pos"] = "32346896"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceDonor)

        # checks that 4th to last base in exon is NOT counted as in splice donor
        self.variant["Pos"] = "32370554"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceDonor)

        # checks that region after  exon 27 is NOT counted as in splice donor
        self.variant["Pos"] = "32399672"
        inSpliceDonor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceDonor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    def test_getVarSpliceRegionBoundsDonorBRCA2(self, varInSpliceRegion, getRefSpliceDonorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position) 
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = True
        deNovo = False

        # checks that variant in exon 16 splice donor region boundaries are returned correctly
        self.variant["Pos"] = "32356608"
        spliceDonorRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor, deNovo=deNovo)
        self.assertEquals(spliceDonorRegion["exonName"], "exon15")
        self.assertEquals(exonDonorBoundsBRCA2["exon15"]["donorStart"], spliceDonorRegion["donorStart"])
        self.assertEquals(exonDonorBoundsBRCA2["exon15"]["donorEnd"], spliceDonorRegion["donorEnd"])

    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_varInSpliceRegionAcceptorBRCA1(self, getSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Variant is correctly identified as in or NOT in a splice acceptor region
           uses multiple positions across multiple exons for BRCA1
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        donor = False
        deNovo = False

        # checks that -21st base in intron is NOT counted as in splice acceptor
        self.variant["Pos"] = "43067716"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceAcceptor)
        
        # checks that first base in intron counted as in splice acceptor
        self.variant["Pos"] = "43124135"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in intron counted as in splice acceptor
        self.variant["Pos"] = "43115787"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in intron counted as in splice acceptor
        self.variant["Pos"] = "43106534"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that first base in exon counted as in splice acceptor
        self.variant["Pos"] = "43104956"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in exon counted as in splice acceptor
        self.variant["Pos"] = "43104260"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in exon counted as in splice acceptor
        self.variant["Pos"] = "43099878"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that 4th base in exon is NOT counted as in splice acceptor
        self.variant["Pos"] = "43063948"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceAcceptor)

        # checks that region before exon 1 is NOT counted as in splice acceptor
        self.variant["Pos"] = "431254483"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceAcceptor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_getVarSpliceRegionBoundsAcceptorBRCA1(self, varInSpliceRegion, getSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position) 
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        donor = False
        deNovo = False

        # checks that variant in exon 21 splice acceptor region boundaries are returned correctly
        self.variant["Pos"] = "43051117"
        spliceAccRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor, deNovo=deNovo)
        self.assertEquals(spliceAccRegion["exonName"], "exon21")
        self.assertEquals(exonAcceptorBoundsBRCA1["exon21"]["acceptorStart"], spliceAccRegion["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA1["exon21"]["acceptorEnd"], spliceAccRegion["acceptorEnd"])

    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_varInSpliceRegionAcceptorBRCA2(self, getSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Variant is correctly identified as in or NOT in a splice acceptor region
           uses multiple positions across multiple exons for BRCA2
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = False
        deNovo = False

        # checks that -21st base in intron is NOT counted as in splice acceptor
        self.variant["Pos"] = "32357721"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceAcceptor)

        # checks that first base in intron counted as in splice acceptor
        self.variant["Pos"] = "32316402"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in intron counted as in splice acceptor
        self.variant["Pos"] = "32319069"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in intron counted as in splice acceptor
        self.variant["Pos"] = "32325075"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that first base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326101"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326243"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326501"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertTrue(inSpliceAcceptor)

        # checks that 4th base in exon is NOT counted as in splice acceptor
        self.variant["Pos"] = "32362526"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceAcceptor)

        # checks that region before exon 1 is NOT counted as in splice acceptor
        self.variant["Pos"] = "32315479"
        inSpliceAcceptor = calcVarPriors.varInSpliceRegion(self.variant, donor=donor, deNovo=deNovo)
        self.assertFalse(inSpliceAcceptor)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_getVarSpliceRegionBoundsAcceptorBRCA2(self, varInSpliceRegion, getSpliceAcceptorBoundaries):
        '''
        Tests that:
        1. Function returns correct donor boundaries for a given variant (genomic position) 
        '''
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        donor = False
        deNovo = False

        # checks that variant in exon 20 splice acceptor region boundaries are returned correctly
        self.variant["Pos"] = "32370948"
        spliceAccRegion = calcVarPriors.getVarSpliceRegionBounds(self.variant, donor=donor, deNovo=deNovo)
        self.assertEquals(spliceAccRegion["exonName"], "exon20")
        self.assertEquals(exonAcceptorBoundsBRCA2["exon20"]["acceptorStart"], spliceAccRegion["acceptorStart"])
        self.assertEquals(exonAcceptorBoundsBRCA2["exon20"]["acceptorEnd"], spliceAccRegion["acceptorEnd"])

    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCIDomainEnigmaBRCA1(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA1 as defined by ENIGMA rules'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks variant in BRCA 1 RING domain is identified as in ENIGMA CI domain
        self.variant["Pos"] = "43124089"
        inEnigmaCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inEnigmaCI)

        # checks variant in BRCA1 BRCT domain is identified as in ENIGMA CI domain
        self.variant["Pos"] = "43070945"
        inEnigmaCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inEnigmaCI)

        # checks variant NOT in BRCA1 CI domain is NOT identified as in ENIGMA CI domain
        self.variant["Pos"] = "43097274"
        inEnigmaCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertFalse(inEnigmaCI)
        
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCIDomainEnigmaBRCA2(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA2 as defined by ENIGMA rules'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks variant in BRCA2 DNB domain is identified as in ENIGMA CI domain
        self.variant["Pos"] = "32379809"
        inEnigmaCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inEnigmaCI)

        # checks variant NOT in BRCA2 CI domain is NOT identified as in ENIGMA CI domain
        self.variant["Pos"] = "32398354"
        inEnigmaCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertFalse(inEnigmaCI)

    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCIDomainPriorsBRCA1(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA1 as defined by PRIORS webiste'''
        boundaries = "priors"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks variant in BRCA1 initiation domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "43124096"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA1 RING domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "43115746"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA1 BRCT domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "43057092"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks that variant NOT in BRCA1 CI domain is NOT identified as in PRIORS CI domain
        self.variant["Pos"] = "43124090"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertFalse(inPriorsCI)
        
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    def test_varInCIDomainPriorsBRCA2(self, getVarStrand, varInExon):
        '''Tests that variant is correctly identified as in or NOT in CI domain in BRCA2 as defined by PRIORS webiste'''
        boundaries = "priors"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks variant in BRCA2 initiation domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32316462"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA2 PALB2 domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32319092"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in BRCA2 DNB domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32362561"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks variant in TR2/RAD5 domain is identified as in PRIORS CI domain
        self.variant["Pos"] = "32398406"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
        self.assertTrue(inPriorsCI)

        # checks that variant NOT in BRCA2 CI domain is NOT identified as in PRIORS CI domain
        self.variant["Pos"] = "32336283"
        inPriorsCI = calcVarPriors.varInCIDomain(self.variant, boundaries)
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
    @mock.patch('calcVarPriors.varInCIDomain', return_value = True)
    def test_getVarLocationCIDomain(self, varOutsideBoundaries, varInExon, varInSpliceRegion, varInCIDomain):
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
    @mock.patch('calcVarPriors.varInCIDomain', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca1RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_getVarLocationCIDomainSpliceRegionBRCA1(self, varOutsideBoundaries, varInExon, varInCIDomain,
                                                     getRefSpliceDonorBoundaries, getSpliceAcceptorBoundaries):
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
    @mock.patch('calcVarPriors.varInCIDomain', return_value = True)
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_getVarLocationCIDomainSpliceRegionBRCA2(self, varOutsideBoundaries, varInExon, varInCIDomain,
                                                     getRefSpliceDonorBoundaries, getSpliceAcceptorBoundaries):
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
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    def test_getVarLocationSpliceRegionBRCA1(self, varOutsideBoundaries, getExonBoundaries, getRefSpliceDonorBoundaries,
                                             getSpliceAcceptorBoundaries):
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
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    def test_getVarLocationSpliceRegionBRCA2(self, varOutsideBoundaries, getExonBoundaries, getRefSpliceDonorBoundaries,
                                             getSpliceAcceptorBoundaries):
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
    @mock.patch('calcVarPriors.varInCIDomain', return_value = False)
    def test_getVarLocationInExon(self, varOutsideBoundaries, varInExon, varInSpliceRegion, varInCIDomain):
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

    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {"inExonicPortion": True})
    def test_varInExonicPortionTrue(self, getMaxMaxEntScanScoreSlidingWindowSNS):
        '''Tests that varInExonicPortion returns True if variant is in exonic portion of window'''
        inExonicPortion = calcVarPriors.varInExonicPortion(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                           donor=True, deNovo=False)
        self.assertTrue(inExonicPortion)

    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {"inExonicPortion": False})
    def test_varInExonicPortionFalse(self, getMaxMaxEntScanScoreSlidingWindowSNS):
        '''Tests that varInExonicPortion returns False if variant is NOT in exonic portion of window'''
        inExonicPortion = calcVarPriors.varInExonicPortion(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH,
                                                           donor=False, deNovo=False)
        self.assertFalse(inExonicPortion)

    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {"varWindowPosition": 2})
    def test_getVarWindowPositionDonorFirstThree(self, getMaxMaxEntScanScoreSlidingWindowSNS):
        '''Tests that function returns correct value for variant in first 3 bp of window'''
        windowPos = calcVarPriors.getVarWindowPosition(self.variant, donor=True, deNovo=False)
        self.assertEquals(windowPos, 2)

    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {"varWindowPosition": 7})
    def test_getVarWindowPositionDonorLastSix(self, getMaxMaxEntScanScoreSlidingWindowSNS):
        '''Tests that function returns correct value for variant after first 3 bp of window'''
        windowPos = calcVarPriors.getVarWindowPosition(self.variant, donor=True, deNovo=False)
        self.assertEquals(windowPos, 7)

    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {"varWindowPosition": 18})
    def test_getVarWindowPositionAcceptor(self, getMaxMaxEntScanScoreSlidingWindowSNS):
        '''Tests that functions returns correct value for de novo acceptor variant'''
        windowPos = calcVarPriors.getVarWindowPosition(self.variant, donor=False, deNovo=True)
        self.assertEquals(windowPos, 18)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon18")
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca1RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "tctgtaagt")
    @mock.patch('calcMaxEntScanMeanStd.runMaxEntScan', return_value = 7.96)
    @mock.patch('calcVarPriors.getZScore', return_value = 0.009407098111019249)
    def test_getClosestSpliceSiteScoresInExonDonorBRCA1(self, varInExon, getVarExonNumberSNS, getRefSpliceDonorBoundaries,
                                                        varInSpliceRegion, getFastaSeq, runMaxEntScan, getZScore):
        '''Tests function for variant in exon to get closest splice donor site in minus strand gene (BRCA1)'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["HGVS_cDNA"] = "c.5104A>G"
        self.variant["Pos"] = "43063922"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        deNovoOffset = 0
        closestScores = calcVarPriors.getClosestSpliceSiteScores(self.variant, deNovoOffset, donor=True, deNovo=False, accDonor=True)
        self.assertEquals(closestScores["exonName"], "exon18")

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon8")
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca2RefSpliceAcceptorBounds)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "CATAAATTTTTATCTTACAGTCA")
    @mock.patch('calcMaxEntScanMeanStd.runMaxEntScan', return_value = 8.03)
    @mock.patch('calcVarPriors.getZScore', return_value = 0.018528005635431572)
    def test_getClosestSpliceSiteScoresInExonAccBRCA2(self, varInExon, getVarExonNumberSNS, getSpliceAcceptorBoundaries,
                                                      varInSpliceRegion, getFastaSeq, runMaxEntScan, getZScore):
        '''Tests function for variant in exon to get closest splice acceptor site in plus strand gene (BRCA2)'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["HGVS_cDNA"] = "c.648A>G"
        self.variant["Pos"] = "32329459"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "G"
        deNovoOffset = 0
        closestScores = calcVarPriors.getClosestSpliceSiteScores(self.variant, deNovoOffset, donor=False, deNovo=False, accDonor=False)
        self.assertEquals(closestScores["exonName"], "exon8")

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon20")
    @mock.patch('calcVarPriors.getRefSpliceDonorBoundaries', return_value = brca2RefSpliceDonorBounds)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = {'exonName': 'exon20',
                                                                          'donorStart': 32371098,
                                                                          'donorEnd': 32371106})
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "AAGGTAAAA")
    @mock.patch('calcMaxEntScanMeanStd.runMaxEntScan', return_value = 8.38)
    @mock.patch('calcVarPriors.getZScore', return_value = 0.18974233990730782)
    def test_getClosestSpliceSiteScoresInRefDonorExonicBRCA2(self, varInExon, getVarExonNumberSNS, getRefSpliceDonorBoundaries,
                                                             varInSpliceRegion, getVarSpliceRegionBounds, getFastaSeq,
                                                             runMaxEntScan, getZScore):
        '''Tests function for variant in exonic portion of ref donor site to get closest splice donor site in plus strand gene (BRCA2)'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["HGVS_cDNA"] = "c.8631A>T"
        self.variant["Pos"] = "32371099"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        deNovoOffset = 0
        closestScores = calcVarPriors.getClosestSpliceSiteScores(self.variant, deNovoOffset, donor=True, deNovo=False, accDonor=False)
        self.assertEquals(closestScores["exonName"], "exon20")

    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = {'acceptorStart': 43106553,
                                                                          'exonName': 'exon5',
                                                                          'acceptorEnd': 43106531})
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "tctttctttataatttatagatt")
    @mock.patch('calcMaxEntScanMeanStd.runMaxEntScan', return_value = 8.19)
    @mock.patch('calcVarPriors.getZScore', return_value = 0.08427254176115757)
    def test_getClosestSpliceSiteScoresInRefAccIntronicBRCA1(self, varInExon, varInSpliceRegion, getVarSpliceRegionBounds,
                                                             getFastaSeq, runMaxEntScan, getZScore):
        '''
        Tests function for variant in intronic portion of ref acceptor site 
        to get closest splice acceptor site in minus strand gene (BRCA1)
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["HGVS_cDNA"] = "c.135-11a>C"
        self.variant["Pos"] = "43106544"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        deNovoOffset = 0
        closestScores = calcVarPriors.getClosestSpliceSiteScores(self.variant, deNovoOffset, donor=False, deNovo=False, accDonor=False)
        self.assertEquals(closestScores["exonName"], "exon5")

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarSpliceRegionBounds', return_value = {'acceptorStart': 32325056,
                                                                          'exonName': 'exon4',
                                                                          'acceptorEnd': 32325085})
    @mock.patch('calcVarPriors.getFastaSeq', return_value = "GAATTATTGTACTGTTTCAGGAA")
    @mock.patch('calcMaxEntScanMeanStd.runMaxEntScan', return_value = 8.41)
    @mock.patch('calcVarPriors.getZScore', return_value = 0.174671278934031)
    def test_getClosestSpliceSiteScoresInDeNovoAccBRCA2(self, varInExon, varInSpliceRegion, getVarSpliceRegionBounds,
                                                        getFastaSeq, runMaxEntScan, getZScore):
        '''Tests function for variant in exon to get closest splice acceptor site for de novo acceptor in plus strand gene (BRCA2)'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["HGVS_cDNA"] = "c.320G>C"
        self.variant["Pos"] = "32325079"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        deNovoOffset = STD_DE_NOVO_OFFSET
        closestScores = calcVarPriors.getClosestSpliceSiteScores(self.variant, deNovoOffset, donor=False, deNovo=True, accDonor=False)
        self.assertEquals(closestScores["exonName"], "exon4")


    def test_isCIDomainInRegionBRCA1(self):
        '''
        Tests that region overlap is identified correctly for a variant on minus strand gene (BRCA1)
        '''
        self.variant["Gene_Symbol"] = "BRCA1"

        boundaries = "enigma"
        # region that includes ENIGMA BRCT domain
        regionStart = 43067625
        regionEnd = 43063950
        CIDomainInRegion = calcVarPriors.isCIDomainInRegion(regionStart, regionEnd, boundaries, self.variant["Gene_Symbol"])
        self.assertTrue(CIDomainInRegion)

        # region that does not include any PRIORS CI domains
        boundaries = "priors"
        regionStart = 43095923
        regionEnd = 43095857
        CIDomainInRegion = calcVarPriors.isCIDomainInRegion(regionStart, regionEnd, boundaries, self.variant["Gene_Symbol"])
        self.assertFalse(CIDomainInRegion)

    def test_isCIDomainInRegionBRCA2(self):
        '''
        Tests that region overlap is identified correctly for a variant on plus strand gene (BRCA2)
        '''
        self.variant["Gene_Symbol"] = "BRCA2"

        # region that does not include any ENIGMA CI domains
        boundaries = "enigma"
        regionStart = 32319089
        regionEnd = 32325063
        CIDomainInRegion = calcVarPriors.isCIDomainInRegion(regionStart, regionEnd, boundaries, self.variant["Gene_Symbol"])
        self.assertFalse(CIDomainInRegion)

        # region that includeds PRIORS DNB domain
        boundaries = "priors"
        regionStart = 32379502
        regionEnd = 32379751
        CIDomainInRegion = calcVarPriors.isCIDomainInRegion(regionStart, regionEnd, boundaries, self.variant["Gene_Symbol"])
        self.assertTrue(CIDomainInRegion)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon9")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    def test_getRefExonLengthBRCA1(self, varInExon, getVarExonNumberSNS, getExonBoundaries, getVarStrand):
        '''Tests that exon length is correctly calculated for minus strand (BRCA1) exon'''
        self.variant["Gene_Symbol"] = "BRCA1"
        exon9PlusSeq = "CTGCAATAAGTTGCCTTATTAACGGTATCTTCAGAAGAATCAGATC"
        refExonLength = calcVarPriors.getRefExonLength(self.variant)
        self.assertEquals(refExonLength, len(exon9PlusSeq))

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon5")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    def test_getRefExonLengthBRCA2(self, varInExon, getVarExonNumberSNS, getExonBoundaries, getVarStrand):
        '''Tests that exon length is correctly calculated for plus strand (BRCA2) exon'''
        self.variant["Gene_Symbol"] = "BRCA2"
        exon5PlusSeq = "TCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAG"
        refExonLength = calcVarPriors.getRefExonLength(self.variant)
        self.assertEquals(refExonLength, len(exon5PlusSeq))

    def test_getNewSplicePositionBRCA1InExonicPortion(self):
        '''Tests that new splice position is calculated correctly for minus strand (BRCA1) variant with max MES in exonic portion'''
        varStrand = "-"
        inExonicPortion = True
        varGenPos = "43104189"
        varWindowPos = 3
        newSplicePos = calcVarPriors.getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, STD_EXONIC_PORTION)
        # because varWindowPos == 3, cut will occur after variant
        actualNewSplicePos = 43104189
        self.assertEquals(newSplicePos, actualNewSplicePos)
        
    def test_getNewSplicePositionBRCA1NotInExonicPortion(self):
        '''Tests that new splice position is calculated correctly for minus strand (BRCA1) variant with max MES NOT in exonic portion'''
        varStrand = "-"
        inExonicPortion = False
        varGenPos = "43104249"
        varWindowPos = 6
        newSplicePos = calcVarPriors.getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, STD_EXONIC_PORTION)
        # because varWindowPos == 6, cut will occur 3 bases to the left of the variant
        actualNewSplicePos = 43104252
        self.assertEquals(newSplicePos, actualNewSplicePos)
        
    def test_getNewSplicePositionBRCA2InExonicPortion(self):
        '''Tests that new splice position is calculated correctly for plus strand (BRCA2) variant with max MES in exonic portion'''
        varStrand = "+"
        inExonicPortion = True
        varGenPos = "32354881"
        varWindowPos = 2
        newSplicePos = calcVarPriors.getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, STD_EXONIC_PORTION)
        # because varWindowPos == 2, cut will occur 1 base to the right of the variant
        actualNewSplicePos = 32354882
        self.assertEquals(newSplicePos, actualNewSplicePos)
        
    def test_getNewSplicePositionBRCA2NotInExonicPortion(self):
        '''Tests that new splice position is calculated correctly for plus strand (BRCA2) variant with max MES NOT in exonic portion'''
        varStrand = "+"
        inExonicPortion = False
        varGenPos = "32326277"
        varWindowPos = 8
        newSplicePos = calcVarPriors.getNewSplicePosition(varGenPos, varStrand, varWindowPos, inExonicPortion, STD_EXONIC_PORTION)
        # because varWindowPos == 8, cut will occur 5 bases to the left of the variant
        actualNewSplicePos = 32326272
        self.assertEquals(newSplicePos, actualNewSplicePos)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon21")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -9.9,
                                                                                       'altMaxEntScanScore': -7.45,
                                                                                       'altZScore': -6.6071787973194605,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 2,
                                                                                       'refZScore': -7.659134374464476})
    @mock.patch('calcVarPriors.getNewSplicePosition', return_value = 43051109)
    def test_getAltExonLengthBRCA1(self, varInExon, getVarExonNumberSNS, getExonBoundaries,
                                   getMaxMaxEntScanScoreSlidingWindowSNS, getNewSplicePosition):
        '''Tests that exon length is correctly calculated for minus strand (BRCA1) exon'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["HGVS_cDNA"] = "c.5285G>C"
        expectedCutSeq = "ATCTTCACG"
        altExonLength = calcVarPriors.getAltExonLength(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(altExonLength, len(expectedCutSeq))

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon13")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -9.89,
                                                                                       'altMaxEntScanScore': -1.7,
                                                                                       'altZScore': -4.13830346320361,
                                                                                       'varWindowPosition': 4,
                                                                                       'inExonicPortion': False,
                                                                                       'refZScore': -7.65484067823123})
    @mock.patch('calcVarPriors.getNewSplicePosition', return_value = 32346831)
    def test_getAltExonLengthBRCA2(self, varInExon, getVarExonNumberSNS, getExonBoundaries,
                                   getMaxMaxEntScanScoreSlidingWindowSNS, getNewSplicePosition):
        '''Tests that exon length is correctly calculated for plus strand (BRCA2) exon'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["HGVS_cDNA"] = "c.6943A>G"
        expectedCutSeq = "GCACA"
        altExonLength = calcVarPriors.getAltExonLength(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(altExonLength, len(expectedCutSeq))

    def test_compareRefAltExonLengths(self):
        '''Tests that function correctly determines if ref and alt exons are in same reading frame'''
        # ref and alt exons that share the same reading frame
        refLength = 45
        altLength = 33
        inFrame = calcVarPriors.compareRefAltExonLengths(refLength, altLength)
        self.assertTrue(inFrame)

        # ref and alt exons that do NOT share the smae reading frame
        refLength = 162
        altLength = 103
        inFrame = calcVarPriors.compareRefAltExonLengths(refLength, altLength)
        self.assertFalse(inFrame)

    @mock.patch('calcVarPriors.getRefExonLength', return_value = 45)
    @mock.patch('calcVarPriors.getAltExonLength', return_value = 30)
    @mock.patch('calcVarPriors.compareRefAltExonLengths', return_value = True)    
    def test_isSplicingWindowInFrameTrue(self, getRefExonLength, getAltExonLength, compareRefAltExonLengths):
        '''Tests that if splicing window is in frame, function returns true'''
        inFrame = calcVarPriors.isSplicingWindowInFrame(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertTrue(inFrame)
        
    @mock.patch('calcVarPriors.getRefExonLength', return_value = 45)
    @mock.patch('calcVarPriors.getAltExonLength', return_value = 29)
    @mock.patch('calcVarPriors.compareRefAltExonLengths', return_value = False)    
    def test_isSplicingWindowInFrameFalse(self, getRefExonLength, getAltExonLength, compareRefAltExonLengths):
        '''Tests that if splicing window is NOT in frame, function returns false'''
        inFrame = calcVarPriors.isSplicingWindowInFrame(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertFalse(inFrame)        

    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon16")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -7.31,
                                                                                       'altMaxEntScanScore': -7.32,
                                                                                       'altZScore': -6.551360746287276,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 1,
                                                                                       'refZScore': -6.547067050054031})
    @mock.patch('calcVarPriors.getNewSplicePosition', return_value = 43070934)
    def test_compareDeNovoWildTypeSplicePosBRCA1True(self, getVarStrand, getVarExonNumberSNS, getExonBoundaries,
                                                     getMaxMaxEntScanScoreSlidingWindowSNS, getNewSplicePosition):
        '''
        Tests that comparsion between de novo and wild-type splice position is correct for a BRCA1 variant that:
            1. has highest scoring window with variant in first three nucleotides
            2. AND distance between de novo and wild-type donor is divisble by 3
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43070936"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "T"
        isDivisible = calcVarPriors.compareDeNovoWildTypeSplicePos(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertTrue(isDivisible)
        
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon9")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca1Exons)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -2.95,
                                                                                       'altMaxEntScanScore': 5.56,
                                                                                       'altZScore': -1.0210799978677707,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 4,
                                                                                       'refZScore': -4.67501549235923})
    @mock.patch('calcVarPriors.getNewSplicePosition', return_value = 43097266)
    def test_compareDeNovoWildTypeSplicePosBRCA1False(self, getVarStrand, getVarExonNumberSNS, getExonBoundaries,
                                                      getMaxMaxEntScanScoreSlidingWindowSNS, getNewSplicePosition):
        '''
        Tests that comparsion between de novo and wild-type splice position is correct for a BRCA1 variant that:
            1. has highest scoring window with variant in last six nucleotides
            2. AND distance between de novo and wild-type donor is NOT divisble by 3
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43097265"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        isDivisible = calcVarPriors.compareDeNovoWildTypeSplicePos(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertFalse(isDivisible)
                
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon4")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -9.22,
                                                                                       'altMaxEntScanScore': -1.57,
                                                                                       'altZScore': -4.082485412171425,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 5,
                                                                                       'refZScore': -7.367163030603819})
    @mock.patch('calcVarPriors.getNewSplicePosition', return_value = 32325178)
    def test_compareDeNovoWildTypeSplicePosBRCA2True(self, getVarStrand, getVarExonNumberSNS, getExonBoundaries,
                                                     getMaxMaxEntScanScoreSlidingWindowSNS, getNewSplicePosition):
        '''
        Tests that comparsion between de novo and wild-type splice position is correct for a BRCA2 variant that:
            1. has highest scoring window with variant in last six nucleotides
            2. AND distance between de novo and wild-type donor is divisble by 3
        '''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32325180"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        isDivisible = calcVarPriors.compareDeNovoWildTypeSplicePos(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertTrue(isDivisible)

    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon23")
    @mock.patch('calcVarPriors.getExonBoundaries', return_value = brca2Exons)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -3.67,
                                                                                       'altMaxEntScanScore': -0.32,
                                                                                       'altZScore': -3.5457733830158054,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 2,
                                                                                       'refZScore': -4.984161621152867})
    @mock.patch('calcVarPriors.getNewSplicePosition', return_value = 32379873)
    def test_compareDeNovoWildTypeSplicePosBRCA2False(self, getVarStrand, getVarExonNumberSNS, getExonBoundaries,
                                                      getMaxMaxEntScanScoreSlidingWindowSNS, getNewSplicePosition):
        '''
        Tests that comparsion between de novo and wild-type splice position is correct for a BRCA2 variant that:
            1. has highest scoring window with variant in first three nucleotides
            2. AND distance between de novo and wild-type donor is NOT divisble by 3
        '''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32379872"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "A"
        isDivisible = calcVarPriors.compareDeNovoWildTypeSplicePos(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertFalse(isDivisible)

    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInExonicPortion', return_value = True)
    def test_getPriorProbSpliceRescueNonsenseSNSInExonicPortion(self, getVarConsequences, varInExon, varInExonicPortion):
        '''Tests that variant in exonic portion of highest scoring window is assigned correct prior prob and splice rescue flag'''
        boundaries = "enigma"
        spliceRescueInfo = calcVarPriors.getPriorProbSpliceRescueNonsenseSNS(self.variant, boundaries, accDonor=False)
        self.assertEquals(spliceRescueInfo["priorProb"], priorProbs["high"])
        self.assertEquals(spliceRescueInfo["enigmaClass"], enigmaClasses["class4"])
        self.assertEquals(spliceRescueInfo["spliceRescue"], 0)
        self.assertEquals(spliceRescueInfo["spliceFlag"], 0)
        self.assertEquals(spliceRescueInfo["frameshift"], 0)

    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInExonicPortion', return_value = False)
    @mock.patch('calcVarPriors.isSplicingWindowInFrame', return_value = False)
    def test_getPriorProbSpliceRescueNonsenseSNSFrameshift(self, getVarConsequences, varInExon, varInExonicPortion, isSplicingWindowInFrame):
        '''Tests that variant that causes a frameshift is assigned correct prior prob and splice rescue flag'''
        boundaries = "enigma"
        spliceRescueInfo = calcVarPriors.getPriorProbSpliceRescueNonsenseSNS(self.variant, boundaries, accDonor=False)
        self.assertEquals(spliceRescueInfo["priorProb"], priorProbs["pathogenic"])
        self.assertEquals(spliceRescueInfo["enigmaClass"], enigmaClasses["class5"])
        self.assertEquals(spliceRescueInfo["spliceRescue"], 0)
        self.assertEquals(spliceRescueInfo["spliceFlag"], 0)
        self.assertEquals(spliceRescueInfo["frameshift"], 1)
        
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInExonicPortion', return_value = False)
    @mock.patch('calcVarPriors.isSplicingWindowInFrame', return_value = True)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon20")
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    @mock.patch('calcVarPriors.getVarWindowPosition', return_value = 6)
    @mock.patch('calcVarPriors.isCIDomainInRegion', return_value = True)
    @mock.patch('calcVarPriors.compareDeNovoWildTypeSplicePos', return_value = True)
    def test_getPriorProbSpliceRescueNonsenseSNSCIRegionEnigma(self, getVarConsequences, varInExon, varInExonicPortion,
                                                               isSplicingWindowInFrame, getVarStrand, getVarExonNumberSNS,
                                                               getSpliceAcceptorBoundaries, getVarWindowPosition,
                                                               isCIDomainInRegion, compareDeNovoWildTypeSplicePos):
        '''Tests that variant that truncates part of ENGIMA CI domain is assigned correct prior prob and splice rescue flag'''
        boundaries = "enigma"
        spliceRescueInfo = calcVarPriors.getPriorProbSpliceRescueNonsenseSNS(self.variant, boundaries, accDonor=False)
        self.assertEquals(spliceRescueInfo["priorProb"], priorProbs["pathogenic"])
        self.assertEquals(spliceRescueInfo["enigmaClass"], enigmaClasses["class5"])
        self.assertEquals(spliceRescueInfo["spliceRescue"], 0)
        self.assertEquals(spliceRescueInfo["spliceFlag"], 0)
        self.assertEquals(spliceRescueInfo["frameshift"], 0)
        
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInExonicPortion', return_value = False)
    @mock.patch('calcVarPriors.isSplicingWindowInFrame', return_value = True)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon2")
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    @mock.patch('calcVarPriors.getVarWindowPosition', return_value = 7)
    @mock.patch('calcVarPriors.isCIDomainInRegion', return_value = True)
    @mock.patch('calcVarPriors.compareDeNovoWildTypeSplicePos', return_value = True)
    def test_getPriorProbSpliceRescueNonsenseSNSCIRegionPriors(self, getVarConsequences, varInExon, varInExonicPortion,
                                                               isSplicingWindowInFrame, getVarStrand, getVarExonNumberSNS,
                                                               getSpliceAcceptorBoundaries, getVarWindowPosition,
                                                               isCIDomainInRegion, compareDeNovoWildTypeSplicePos):
        '''Tests that variant that truncates part of PRIORS CI domain is assigned correct prior prob and splice rescue flag'''
        boundaries = "priors"
        spliceRescueInfo = calcVarPriors.getPriorProbSpliceRescueNonsenseSNS(self.variant, boundaries, accDonor=False)
        self.assertEquals(spliceRescueInfo["priorProb"], priorProbs["pathogenic"])
        self.assertEquals(spliceRescueInfo["enigmaClass"], enigmaClasses["class5"])
        self.assertEquals(spliceRescueInfo["spliceRescue"], 0)
        self.assertEquals(spliceRescueInfo["spliceFlag"], 0)
        self.assertEquals(spliceRescueInfo["frameshift"], 0)
        
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInExonicPortion', return_value = False)
    @mock.patch('calcVarPriors.isSplicingWindowInFrame', return_value = True)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "+")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon6")
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    @mock.patch('calcVarPriors.getVarWindowPosition', return_value = 4)
    @mock.patch('calcVarPriors.isCIDomainInRegion', return_value = False)
    @mock.patch('calcVarPriors.compareDeNovoWildTypeSplicePos', return_value = False)
    def test_getPriorProbSpliceRescueNonsenseSNSNotDivisible(self, getVarConsequences, varInExon, varInExonicPortion,
                                                             isSplicingWindowInFrame, getVarStrand, getVarExonNumberSNS,
                                                             getSpliceAcceptorBoundaries, getVarWindowPosition,
                                                             isCIDomainInRegion, compareDeNovoWildTypeSplicePos):
        '''
        Tests that variant that causes a frameshift (due to difference de novo vs wild-type splice position) 
        is assigned correct prior prob and splice rescue flag
        '''
        boundaries = "enigma"
        spliceRescueInfo = calcVarPriors.getPriorProbSpliceRescueNonsenseSNS(self.variant, boundaries, accDonor=False)
        self.assertEquals(spliceRescueInfo["priorProb"], priorProbs["pathogenic"])
        self.assertEquals(spliceRescueInfo["enigmaClass"], enigmaClasses["class5"])
        self.assertEquals(spliceRescueInfo["spliceRescue"], 0)
        self.assertEquals(spliceRescueInfo["spliceFlag"], 0)
        self.assertEquals(spliceRescueInfo["frameshift"], 1)

    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInExonicPortion', return_value = False)
    @mock.patch('calcVarPriors.isSplicingWindowInFrame', return_value = True)
    @mock.patch('calcVarPriors.getVarStrand', return_value = "-")
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon17")
    @mock.patch('calcVarPriors.getSpliceAcceptorBoundaries', return_value = brca1RefSpliceAcceptorBounds)
    @mock.patch('calcVarPriors.getVarWindowPosition', return_value = 5)
    @mock.patch('calcVarPriors.isCIDomainInRegion', return_value = False)
    @mock.patch('calcVarPriors.compareDeNovoWildTypeSplicePos', return_value = True)
    def test_getPriorProbSpliceRescueNonsenseSNSWithSpliceFlag(self, getVarConsequences, varInExon, varInExonicPortion,
                                                               isSplicingWindowInFrame, getVarStrand, getVarExonNumberSNS,
                                                               getSpliceAcceptorBoundaries, getVarWindowPosition,
                                                               isCIDomainInRegion, compareDeNovoWildTypeSplicePos):
        '''Tests that variant with possibility of splice rescue is assigned correct splice rescue and splicing flag'''
        boundaries = "enigma"
        spliceRescueInfo = calcVarPriors.getPriorProbSpliceRescueNonsenseSNS(self.variant, boundaries, accDonor=False)
        self.assertEquals(spliceRescueInfo["priorProb"], priorProbs["NA"])
        self.assertEquals(spliceRescueInfo["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(spliceRescueInfo["spliceRescue"], 1)
        self.assertEquals(spliceRescueInfo["spliceFlag"], 1)
        self.assertEquals(spliceRescueInfo["frameshift"], 0)
        
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
    def test_getPriorProbRefSpliceDonorSNSLowProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                       getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variant that creates a resaonble splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that creates a reasonable splice donor site
        self.variant["Pos"] = "43076489"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSModerateProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variant that weakens a reasonably strong splice donor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that weakens a reasonably strong splice donor site
        self.variant["Pos"] = "43076485"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSHighProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                        getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA1 variant that further weakens a weak splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that further weakens a weak splice donor site
        self.variant["Pos"] = "43104120"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSImprovedProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variant that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that makes a splice donor site stronger or equally strong
        self.variant["Pos"] = "43076490"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSLowProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                       getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that creates a resaonble splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that creates a reasonable splice donor site
        self.variant["Pos"] = "32346895"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSModerateProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that weakens a reasonably strong splice donor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that weakens a reasonably strong splice donor site
        self.variant["Pos"] = "32346899"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSHighProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                        getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA2 variant that further weakens a weak splice donor site'''  
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that further weakens a weak splice donor site
        self.variant["Pos"] = "32362693"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceDonorSNSImprovedProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                            getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that makes a splice donor site stronger or equally strong
        self.variant["Pos"] = "32346902"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbRefSpliceDonorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSLowProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                          getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variants that creates a resaonble splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that creates a reasonable splice acceptor site
        self.variant["Pos"] = "43104275"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSModerateProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                               getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variants that weakens a reasonably strong splice acceptor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that weakens a reasonably strong splice acceptor site
        self.variant["Pos"] = "43104267"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSHighProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                           getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA1 variant that further weakens a weak splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"

        # checks prior prob for BRCA1 variant that further weakens a weak splice acceptor site
        self.variant["Pos"] = "43124116"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSImprovedProbBRCA1(self, getFastaSeq, getVarType, getVarLocation,
                                                               getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA1 variants that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
    
        # checks prior prob for BRCA1 variant that makes a splice acceptor site stronger or equally strong
        self.variant["Pos"] = "43104261"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSLowProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                          getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that creates a resaonble splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that creates a reasonable splice acceptor site
        self.variant["Pos"] = "32379751"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSModerateProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                               getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that weakens a reasonably strong splice acceptor'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that weakens a reasonably strong splice acceptor site
        self.variant["Pos"] = "32379746"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSHighProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                           getVarSpliceRegionBounds, getRefAltScores):
        '''Tests fucntion for BRCA2 variant that further weakens a weak splice acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that further weakens a weak splice acceptor site
        self.variant["Pos"] = "32356427"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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
    def test_getPriorProbRefSpliceAcceptorSNSImprovedProbBRCA2(self, getFastaSeq, getVarType, getVarLocation,
                                                                getVarSpliceRegionBounds, getRefAltScores):
        '''Tests function for BRCA2 variant that makes a splice site stronger or equally strong'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"

        # checks prior prob for BRCA2 variant that makes a splice acceptor site stronger or equally strong
        self.variant["Pos"] = "32379731"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbRefSpliceAcceptorSNS(self.variant, boundaries)
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

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon24")
    def test_varInIneligibleDeNovoExonDonorBRCA1True(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in last exon of BRCA1 is correctly identified as ineligible for de novo donor'''
        self.variant["Gene_Symbol"] = "BRCA1"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=True)
        self.assertTrue(ineligibleDeNovo)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon5")
    def test_varInIneligibleDeNovoExonDonorBRCA1False(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in other exon of BRCA1 is correctly identified as eligible for de novo donor'''
        self.variant["Gene_Symbol"] = "BRCA1"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=True)
        self.assertFalse(ineligibleDeNovo)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon27")
    def test_varInIneligibleDeNovoExonDonorBRCA2True(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in last exon of BRCA2 is correctly identified as ineligible for de novo donor'''
        self.variant["Gene_Symbol"] = "BRCA2"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=True)
        self.assertTrue(ineligibleDeNovo)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon19")
    def test_varInIneligibleDeNovoExonDonorBRCA2False(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in other exon of BRCA2 is correctly identified as eligible for de novo donor'''
        self.variant["Gene_Symbol"] = "BRCA2"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=True)
        self.assertFalse(ineligibleDeNovo)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon1")
    def test_varInIneligibleDeNovoExonAccBRCA1True(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in first exon of BRCA1 is correctly identified as ineligible for de novo acceptor'''
        self.variant["Gene_Symbol"] = "BRCA1"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=False)
        self.assertTrue(ineligibleDeNovo)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon16")
    def test_varInIneligibleDeNovoExonAccBRCA1False(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in other exon of BRCA1 is correctly identified as eligible for de novo acceptor'''
        self.variant["Gene_Symbol"] = "BRCA1"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=False)
        self.assertFalse(ineligibleDeNovo)
        
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon1")
    def test_varInIneligibleDeNovoExonAccBRCA2True(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in first exon of BRCA2 is correctly identified as ineligible for de novo acceptor'''
        self.variant["Gene_Symbol"] = "BRCA2"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=False)
        self.assertTrue(ineligibleDeNovo)

    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getVarExonNumberSNS', return_value = "exon2")
    def test_varInIneligibleDeNovoExonAccBRCA2False(self, varInExon, getVarExonNumberSNS):
        '''Tests that variant in other exon of BRCA2 is correctly identified as eligible for de novo acceptor'''
        self.variant["Gene_Symbol"] = "BRCA2"
        ineligibleDeNovo = calcVarPriors.varInIneligibleDeNovoExon(self.variant, donor=False)
        self.assertFalse(ineligibleDeNovo)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = True)
    def test_getPriorProbDeNovoDonorSNSELastExonBRCA1(self, getVarType, varInExon, varInSpliceRegion,
                                                      varInIneligibleDeNovoExon):
        '''Tests that variant in last exon of BRCA1 is correctly assigned de novo prob'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43045765"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = True)
    def test_getPriorProbDeNovoDonorSNSELastExonBRCA2(self, getVarType, varInExon, varInSpliceRegion,
                                                      varInIneligibleDeNovoExon):
        '''Tests that variant in last exon of BRCA2 is correclty assigned de novo prob'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32398180"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 1.49,
                                                                                       'altMaxEntScanScore': -3.56,
                                                                                       'altZScore': -4.936930962587172,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 2,
                                                                                       'refZScore': -2.7686143647984682})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 1.1729987773204027,
                                                                            'maxEntScanScore': 10.67})
    def test_getPriorProbDeNovoDonorSNSExonRefGreaterAltBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                              getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in exon where ref zscore is greater than alt zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43097273"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    # the below variant is in exonic portion of splice site, just making sure that function works if varInExon == False
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -12.18,
                                                                                       'altMaxEntScanScore': -13.41,
                                                                                       'altZScore': -9.166221752333454,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 1,
                                                                                       'refZScore': -8.638097115644324})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.27990996080545155,
                                                                            'maxEntScanScore': 8.59})
    def test_getPriorProbDeNovoDonorSNSSpliceSiteRefGreaterAltBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                    getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in splice site where ref zscore is greater than alt zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43090945"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -3.85,
                                                                                       'altMaxEntScanScore': -7.01,
                                                                                       'altZScore': -6.418256163056683,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 6,
                                                                                       'refZScore': -5.061448153351276})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -1.3516946078276324,
                                                                            'maxEntScanScore': 4.79})
    def test_getPriorProbDeNovoDonorSNSExonRefGreaterAltBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                              getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in exon where ref zscore is greater than alt zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32344578"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    # the below variant is in exonic portion of splice site, just making sure that function works if varInExon == False
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -14.69,
                                                                                       'altMaxEntScanScore': -17.11,
                                                                                       'altZScore': -10.75488935863409,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 1,
                                                                                       'refZScore': -9.71581487018881})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 1.1128870300549731,
                                                                            'maxEntScanScore': 10.53})
    def test_getPriorProbDeNovoDonorSNSSpliceSiteRefGreaterAltBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                    getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in splice site where ref zscore is greater than alt zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32346895"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -5.67,
                                                                                       'altMaxEntScanScore': -3.44,
                                                                                       'altZScore': -4.885406607788233,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 7,
                                                                                       'refZScore': -5.842900867801858})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 1.1300618149879533,
                                                                            'maxEntScanScore': 10.57})
    def test_getPriorProbDeNovoDonorSNSExonLowProbBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                        getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in exon with expected low (0.02) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43076557"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoLow"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    # the below variant is in exonic portion of splice site, just making sure that function works if varInExon == False
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -4.41,
                                                                                       'altMaxEntScanScore': -2.12,
                                                                                       'altZScore': -4.318638704999898,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 7,
                                                                                       'refZScore': -5.301895142412993})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 1.1729987773204027,
                                                                            'maxEntScanScore': 10.67})
    def test_getPriorProbDeNovoDonorSNSSpliceSiteLowProbBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                              getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in splice site with expected low (0.02) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43097244"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoLow"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -12.7,
                                                                                       'altMaxEntScanScore': -5.05,
                                                                                       'altZScore': -5.57669170134067,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 8,
                                                                                       'refZScore': -8.861369319773063})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.6534615330977633,
                                                                            'maxEntScanScore': 9.46})
    def test_getPriorProbDeNovoDonorSNSExonLowProbBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                        getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in exon with expected low (0.02) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32326110"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoLow"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    # the below variant is in exonic portion of splice site, just making sure that function works if varInExon == False
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -6.24,
                                                                                       'altMaxEntScanScore': -4.4,
                                                                                       'altZScore': -5.297601446179749,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 9,
                                                                                       'refZScore': -6.087641553096821})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 1.164411384853913,
                                                                            'maxEntScanScore': 10.65})
    def test_getPriorProbDeNovoDonorSNSSpliceSiteLowProbBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                              getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in splice site with expected low (0.02) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32397044"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoLow"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class2"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -2.47,
                                                                                       'altMaxEntScanScore': 3.73,
                                                                                       'altZScore': -1.8068264085515982,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 8,
                                                                                       'refZScore': -4.468918073163471})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.9196706995589503,
                                                                            'maxEntScanScore': 10.08})
    def test_getPriorProbDeNovoDonorSNSExonModProbBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                        getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in exon with expected moderate (0.3) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43115740"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    # the below variant is in exonic portion of splice site, just making sure that function works if varInExon == False
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -0.35,
                                                                                       'altMaxEntScanScore': 7.4,
                                                                                       'altZScore': -0.2310398909506982,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 5,
                                                                                       'refZScore': -3.558654471715541})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.9196706995589503,
                                                                            'maxEntScanScore': 10.08})
    def test_getPriorProbDeNovoDonorSNSSpliceSiteModProbBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                              getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in splice site with expected moderate (0.3) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43115728"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -3.5,
                                                                                       'altMaxEntScanScore': 5.0,
                                                                                       'altZScore': -1.2615269869294883,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 4,
                                                                                       'refZScore': -4.911168785187702})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.09528102277591848,
                                                                            'maxEntScanScore': 8.16})
    def test_getPriorProbDeNovoDonorSNSExonModProbBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                        getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in exon with expected moderate (0.3) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32332406"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    # the below variant is in exonic portion of splice site, just making sure that function works if varInExon == False
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -0.91,
                                                                                       'altMaxEntScanScore': 6.84,
                                                                                       'altZScore': -0.47148688001241607,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 5,
                                                                                       'refZScore': -3.799101460777258})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.17686125120757246,
                                                                            'maxEntScanScore': 8.35})
    def test_getPriorProbDeNovoDonorSNSSpliceSiteModProbBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                              getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in splice site with expected moderate (0.3) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32316525"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 2.78,
                                                                                       'altMaxEntScanScore': 8.94,
                                                                                       'altZScore': 0.4301893289690249,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 3,
                                                                                       'refZScore': -2.2147275507098687})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.49030107623445457,
                                                                            'maxEntScanScore': 9.08})
    def test_getPriorProbDeNovoDonorSNSExonHighProbBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                         getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA1 variant in exon with expected high (0.64) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43099839"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 2.14,
                                                                                       'altMaxEntScanScore': 9.79,
                                                                                       'altZScore': 0.7951535087948461,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 5,
                                                                                       'refZScore': -2.4895241096375464})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 1.2545790057520567,
                                                                            'maxEntScanScore': 10.86})
    def test_getPriorProbDeNovoDonorSNSExonHighProbBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                         getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''Tests BRCA2 variant in exon with expected high (0.64) prior prob where alt zscore > ref zscore'''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32379455"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 3.77,
                                                                                       'altMaxEntScanScore': 3.37,
                                                                                       'altZScore': -1.9613994729484163,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 3,
                                                                                       'refZScore': -1.7896516236186177})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -2.021511220213846,
                                                                            'maxEntScanScore': 3.23})
    def test_getPriorProbDeNovoDonorSNSExonLowProbGreaterSubBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                  getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''
        Tests BRCA1 variant in exon with expected low (0.02) prior prob that is promoted to moderate prior prob 
        because alt zscore > subsequent z score
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43104184"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 0.21,
                                                                                       'altMaxEntScanScore': 3.25,
                                                                                       'altZScore': -2.012923827747356,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 7,
                                                                                       'refZScore': -3.318207482653823})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -2.07732927124603,
                                                                            'maxEntScanScore': 3.1})
    def test_getPriorProbDeNovoDonorSNSExonLowProbGreaterSubBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                  getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''
        Tests BRCA2 variant in exon with expected low (0.02) prior prob that is promoted to moderate prior prob 
        because alt zscore > subsequent z score
        '''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32362627"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -0.17,
                                                                                       'altMaxEntScanScore': 7.59,
                                                                                       'altZScore': -0.14945966251904425,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 5,
                                                                                       'refZScore': -3.4813679395171317})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -0.9867304280018111,
                                                                            'maxEntScanScore': 5.64})
    def test_getPriorProbDeNovoDonorSNSExonModProbGreaterSubBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                  getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''
        Tests BRCA1 variant in exon with expected moderate (0.3) prior prob that is promoted to high prior prob 
        because alt zscore > subsequent z score
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Referenec_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43094762"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 6.34,
                                                                                       'altMaxEntScanScore': 6.92,
                                                                                       'altZScore': -0.43713731014645635,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 7,
                                                                                       'refZScore': -0.686171691674664})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -2.07732927124603,
                                                                            'maxEntScanScore': 3.1})
    def test_getPriorProbDeNovoDonorSNSExonModProbGreaterSubBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                  getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''
        Tests BRCA2 variant in exon with expected moderate (0.3) prior prob that is promoted to high prior prob 
        because alt zscore > subsequent z score
        '''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Referenec_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32362546"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 3.77,
                                                                                       'altMaxEntScanScore': 9.26,
                                                                                       'altZScore': 0.5675876084328637,
                                                                                       'inExonicPortion': True,
                                                                                       'varWindowPosition': 3,
                                                                                       'refZScore': -1.7896516236186177})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -2.021511220213846,
                                                                            'maxEntScanScore': 3.23})
    def test_getPriorProbDeNovoDonorSNSExonHighProbGreaterSubBRCA1(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                   getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''
        Tests BRCA1 variant in exon with expected high (0.64) prior prob and alt zscore > subsequent z score
        '''
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43104184"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = False)
    @mock.patch('calcVarPriors.varInIneligibleDeNovoExon', return_value = False)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 0.81,
                                                                                       'altMaxEntScanScore': 8.99,
                                                                                       'altZScore': 0.45165781013525,
                                                                                       'inExonicPortion': False,
                                                                                       'varWindowPosition': 4,
                                                                                       'refZScore': -3.0605857086591257})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.4044271515695557,
                                                                            'maxEntScanScore': 8.88})
    def test_getPriorProbDeNovoDonorSNSExonHighProbGreaterSubBRCA2(self, getVarType, varInExon, varInSpliceRegion, varInIneligibleDeNovoExon,
                                                                   getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        '''
        Tests BRCA2 variant in exon with expected high (0.64) prior prob and alt zscore > subsequent z score
        '''
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32363225"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoDonorSNS(self.variant, STD_EXONIC_PORTION, accDonor=False)
        self.assertEquals(priorProb["priorProb"], priorProbs["deNovoHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])
            
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -3.96,
                                                                                       'altMaxEntScanScore': -4.11,
                                                                                       'altZScore': -4.969838672904024,
                                                                                       'inExonicPortion': 'N/A',
                                                                                       'varWindowPosition': 20,
                                                                                       'refZScore': -4.908203170286155})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.2815061501383356,
                                                                            'maxEntScanScore': 8.67})
    def test_getPriorProbDeNovoAccSNSFalseAltLessRefBRCA1(self, getVarType, varInSpliceRegion,
                                                          getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43049200"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbDeNovoAcceptorSNS(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoAccFlag"], 0)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -11.02,
                                                                                       'altMaxEntScanScore': -11.86,
                                                                                       'altZScore': -8.154339641493873,
                                                                                       'inExonicPortion': 'N/A',
                                                                                       'varWindowPosition': 16,
                                                                                       'refZScore': -7.8091808268338125})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.018528005635431572,
                                                                            'maxEntScanScore': 8.03})
    def test_getPriorProbDeNovoAccSNSFalseAltLessRefBRCA2(self, getVarType, varInSpliceRegion,
                                                          getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32329450"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoAcceptorSNS(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoAccFlag"], 0)

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -0.03,
                                                                                       'altMaxEntScanScore': 0.26,
                                                                                       'altZScore': -3.174191029970134,
                                                                                       'inExonicPortion': 'N/A',
                                                                                       'varWindowPosition': 23,
                                                                                       'refZScore': -3.2933530016980126})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': 0.8280076066834324,
                                                                            'maxEntScanScore': 10.0})
    def test_getPriorProbDeNovoAccSNSFlagAltGreaterRefBRCA1(self, getVarType, varInSpliceRegion,
                                                             getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43095922"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoAcceptorSNS(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoAccFlag"], 1)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': -7.9,
                                                                                       'altMaxEntScanScore': 0.7,
                                                                                       'altZScore': -2.9933935556243876,
                                                                                       'inExonicPortion': 'N/A',
                                                                                       'varWindowPosition': 20,
                                                                                       'refZScore': -6.527162372382157})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -1.5305776268269857,
                                                                            'maxEntScanScore': 4.26})
    def test_getPriorProbDeNovoAccSNSFlagAltGreaterRefBRCA2(self, getVarType, varInSpliceRegion,
                                                             getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32370393"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbDeNovoAcceptorSNS(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoAccFlag"], 1)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 3.17,
                                                                                       'altMaxEntScanScore': 4.73,
                                                                                       'altZScore': -1.3374530519576655,
                                                                                       'inExonicPortion': 'N/A',
                                                                                       'varWindowPosition': 22,
                                                                                       'refZScore': -1.9784622791834936})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -1.4031975880833913,
                                                                            'maxEntScanScore': 4.57})
    def test_getPriorProbDeNovoAccSNSFlagBRCA1(self, getVarType, varInSpliceRegion,
                                               getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43082571"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbDeNovoAcceptorSNS(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoAccFlag"], 1)
        
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getMaxMaxEntScanScoreSlidingWindowSNS', return_value = {'refMaxEntScanScore': 2.86,
                                                                                       'altMaxEntScanScore': 5.26,
                                                                                       'altZScore': -1.1196742760411986,
                                                                                       'inExonicPortion': 'N/A',
                                                                                       'varWindowPosition': 13,
                                                                                       'refZScore': -2.1058423179270873})
    @mock.patch('calcVarPriors.getClosestSpliceSiteScores', return_value = {'zScore': -1.5305776268269857,
                                                                            'maxEntScanScore': 4.26})
    def test_getPriorProbDeNovoAccSNSFlagBRCA2(self, getVarType, varInSpliceRegion,
                                               getMaxMaxEntScanScoreSlidingWindowSNS, getClosestSpliceSiteScores):
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32370408"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbDeNovoAcceptorSNS(self.variant, STD_EXONIC_PORTION, STD_DE_NOVO_LENGTH)
        self.assertEquals(priorProb["priorProb"], priorProbs["NA"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["NA"])
        self.assertEquals(priorProb["deNovoAccFlag"], 1)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getPriorProbRefSpliceDonorSNS', return_value = {'refMaxEntScanScore': 10.08,
                                                                               'altMaxEntScanScore': 8.31,
                                                                               'enigmaClass': 'class_2',
                                                                               'altZScore': 0.159686466274593,
                                                                               'priorProb': 0.04,
                                                                               'spliceSite': 1,
                                                                               'refZScore': 0.9196706995589503})
    @mock.patch('calcVarPriors.getPriorProbDeNovoDonorSNS', return_value = {'refMaxEntScanScore': -0.35,
                                                                            'altMaxEntScanScore': -5.25,
                                                                            'enigmaClass': 'N/A',
                                                                            'altZScore': -5.66256562600557,
                                                                            'priorProb': 'N/A',
                                                                            'deNovoDonorFlag': 0,
                                                                            'refZScore': -3.558654471715541})
    @mock.patch('calcVarPriors.getPriorProbProteinSNS', return_value = {'enigmaClass': 'class_3',
                                                                        'priorProb': 0.29})
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "missense_variant")
    def test_getPriorProbSpliceDonorSNSNoDeNovoBRCA1(self, varInSpliceRegion, getVarType, varInExon,
                                                     getPriorProbRefSpliceDonorSNS, getPriorProbDeNovoDonorSNS,
                                                     getPriorProbProteinSNS, getVarConsequences):
        '''Tests that applicable prior for a variant in a reference splice site is assigned correctly (no de novo splicing)'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Chr"] = "17"
        self.variant["HGVS_cDNA"] = "c.134A>T"
        self.variant["Pos"] = self.variant["Hg38_Start"] = self.variant["Hg38_End"] = "43115726"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries, variantData)
        # checks that variant splice site flag is assigned correctly
        self.assertEquals(priorProb["spliceSite"], 1)
        # checks that variant is NOT flagged as a de novo splice donor or acceptor
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)
        self.assertEquals(priorProb["deNovoAccFlag"], 0)
        # checks that prior prob and enigma class are appropriate based on applicable prior
        self.assertEquals(priorProb["applicablePrior"], priorProbs["proteinMod"])
        self.assertEquals(priorProb["applicableEnigmaClass"], enigmaClasses["class3"])
        # checks that protein prior prob, ref prior prob, and de novo prior prob are set appropriately
        self.assertEquals(priorProb["proteinPrior"], priorProbs["proteinMod"])
        self.assertEquals(priorProb["refDonorPrior"], priorProbs["low"])
        self.assertEquals(priorProb["deNovoDonorPrior"], "N/A")
        # checks that a score is present for reference and de novo donor values
        self.assertNotEquals(priorProb["altRefDonorZ"], "-")
        self.assertNotEquals(priorProb["refDeNovoDonorMES"], "-")
        # checks that scores are not present for ref splice acceptor site or de novo splice acceptor sites
        self.assertEquals(priorProb["altRefAccZ"], "-")
        self.assertEquals(priorProb["refDeNovoAccMES"], "-")
        # checks that splice rescue, splice flag, and frameshift are all equal to zero
        self.assertEquals(priorProb["spliceRescue"], 0)
        self.assertEquals(priorProb["spliceFlag"], 0)
        self.assertEquals(priorProb["frameshift"], 0)
        
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getPriorProbRefSpliceDonorSNS', return_value = {'refMaxEntScanScore': 8.88,
                                                                               'altMaxEntScanScore': 5.83,
                                                                               'enigmaClass': 'class_3',
                                                                               'altZScore': -0.9051501995701567,
                                                                               'priorProb': 0.34,
                                                                               'spliceSite': 1,
                                                                               'refZScore': 0.4044271515695557})
    @mock.patch('calcVarPriors.getPriorProbDeNovoDonorSNS', return_value = {'refMaxEntScanScore': -9.44,
                                                                            'altMaxEntScanScore': 1.56,
                                                                            'enigmaClass': 'class_2',
                                                                            'altZScore': -2.7385584911657537,
                                                                            'priorProb': 0.02,
                                                                            'deNovoDonorFlag': 1,
                                                                            'refZScore': -7.461624347735207})
    @mock.patch('calcVarPriors.getPriorProbProteinSNS', return_value = {'enigmaClass': 'class_2',
                                                                        'priorProb': 0.02})
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "missense_variant")
    def test_getPriorProbSpliceDonorSNSWithDeNovoBRCA2(self, varInSpliceRegion, getVarType, varInExon,
                                                       getPriorProbRefSpliceDonorSNS, getPriorProbDeNovoDonorSNS,
                                                       getPriorProbProteinSNS, getVarConsequences):
        '''Tests that applicable prior for a variant in a reference splice site is assigned correctly (with de novo splicing)'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Chr"] = "13"
        self.variant["HGVS_cDNA"] = "c.515A>G"
        self.variant["Pos"] = self.variant["Hg38_Start"] = self.variant["Hg38_End"] = "32326281"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries, variantData)
        # checks that variant splice site flag and de novo splice flag are assigned correctly
        self.assertEquals(priorProb["spliceSite"], 1)
        self.assertEquals(priorProb["deNovoDonorFlag"], 1)
        # checks that variant is not flagged as a de novo splice acceptor
        self.assertEquals(priorProb["deNovoAccFlag"], 0)
        # checks that prior prob and enigma class are appropriate based on applicable prior
        self.assertEquals(priorProb["applicablePrior"], priorProbs["moderate"])
        self.assertEquals(priorProb["applicableEnigmaClass"], enigmaClasses["class3"])
        # checks that protein prior prob, ref prior prob, and de novo prior prob are set appropriately
        self.assertEquals(priorProb["proteinPrior"], priorProbs["deNovoLow"])
        self.assertEquals(priorProb["refDonorPrior"], priorProbs["moderate"])
        self.assertEquals(priorProb["deNovoDonorPrior"], priorProbs["deNovoLow"])
        # checks that a score is present for reference donor score and de novo splice donor score
        self.assertNotEquals(priorProb["altDeNovoDonorZ"], "-")
        self.assertNotEquals(priorProb["refRefDonorMES"], "-")
        # checks that scores are not present for ref splice acceptor site or de novo splice acceptor sites
        self.assertEquals(priorProb["altDeNovoAccZ"], "-")
        self.assertEquals(priorProb["refRefAccMES"], "-")
        # checks that splice rescue, splice flag, and frameshift are all equal to zero
        self.assertEquals(priorProb["spliceRescue"], 0)
        self.assertEquals(priorProb["spliceFlag"], 0)
        self.assertEquals(priorProb["frameshift"], 0)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getPriorProbRefSpliceDonorSNS', return_value = {'refMaxEntScanScore': 9.66,
                                                                               'altMaxEntScanScore': 6.69,
                                                                               'enigmaClass': 'class_3',
                                                                               'altZScore': -0.5358923235110902,
                                                                               'priorProb': 0.34,
                                                                               'spliceSite': 1,
                                                                               'refZScore': 0.7393354577626622})
    @mock.patch('calcVarPriors.getPriorProbDeNovoDonorSNS', return_value = {'refMaxEntScanScore': -1.72,
                                                                            'altMaxEntScanScore': -0.67,
                                                                            'enigmaClass': 'class_2',
                                                                            'altZScore': -3.6960527511793795,
                                                                            'priorProb': 0.02,
                                                                            'deNovoDonorFlag': 1,
                                                                            'refZScore': -4.1468908556701})
    @mock.patch('calcVarPriors.getPriorProbProteinSNS', return_value = {'enigmaClass': 'class_5',
                                                                        'priorProb': 0.99})
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.getPriorProbSpliceRescueNonsenseSNS', return_value = {'spliceFlag': 0,
                                                                                     'enigmaClass': 'class_5',
                                                                                     'frameshift': 1,
                                                                                     'priorProb': 0.99,
                                                                                     'spliceRescue': 0})
    def test_getPriorProbSpliceDonorSNSNonsenseBRCA2(self, varInSpliceRegion, getVarType, varInExon,
                                                     getPriorProbRefSpliceDonorSNS, getPriorProbDeNovoDonorSNS,
                                                     getPriorProbProteinSNS, getVarConsequences,
                                                     getPriorProbSpliceRescueNonsenseSNS):
        '''Tests function for nonsense variant in reference splice donor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Chr"] = "13"
        self.variant["HGVS_cDNA"] = "c.316G>T"
        self.variant["Pos"] = self.variant["Hg38_Start"] = self.variant["Hg38_End"] = "32319325"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbSpliceDonorSNS(self.variant, boundaries, variantData)
        # checks that variant splice site flag and de novo splice flag are assigned correctly
        self.assertEquals(priorProb["spliceSite"], 1)
        self.assertEquals(priorProb["deNovoDonorFlag"], 1)
        # checks that variant is not flagged as a de novo splice acceptor
        self.assertEquals(priorProb["deNovoAccFlag"], 0)
        # checks that prior prob and enigma class are appropriate based on applicable prior
        self.assertEquals(priorProb["applicablePrior"], priorProbs["pathogenic"])
        self.assertEquals(priorProb["applicableEnigmaClass"], enigmaClasses["class5"])
        # checks that protein prior prob, ref prior prob, and de novo prior prob are set appropriately
        self.assertEquals(priorProb["proteinPrior"], priorProbs["pathogenic"])
        self.assertEquals(priorProb["refDonorPrior"], priorProbs["moderate"])
        self.assertEquals(priorProb["deNovoDonorPrior"], priorProbs["deNovoLow"])
        # checks that a score is present for reference donor score and de novo splice donor score
        self.assertNotEquals(priorProb["altDeNovoDonorZ"], "-")
        self.assertNotEquals(priorProb["refRefDonorMES"], "-")
        # checks that scores are not present for ref splice acceptor site or de novo splice acceptor sites
        self.assertEquals(priorProb["altDeNovoAccZ"], "-")
        self.assertEquals(priorProb["refRefAccMES"], "-")
        # checks that splice rescue and splice flag are equal to zero
        self.assertEquals(priorProb["spliceRescue"], 0)
        self.assertEquals(priorProb["spliceFlag"], 0)
        # checks that frameshift is equal to 1, because nonsense variant causes frameshift mutation
        self.assertEquals(priorProb["frameshift"], 1)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getPriorProbRefSpliceAcceptorSNS', return_value = {'refMaxEntScanScore': 10.37,
                                                                                  'altMaxEntScanScore': 11.62,
                                                                                  'enigmaClass': 'class_2',
                                                                                  'altZScore': 1.4936710349564073,
                                                                                  'priorProb': 0.04,
                                                                                  'spliceSite': 1,
                                                                                  'refZScore': 0.9800418464741734})
    @mock.patch('calcVarPriors.getPriorProbDeNovoAcceptorSNS', return_value = {'refMaxEntScanScore': -13.96,
                                                                               'altMaxEntScanScore': -5.9,
                                                                               'enigmaClass': 'N/A',
                                                                               'altZScore': -5.705355670810582,
                                                                               'priorProb': 'N/A',
                                                                               'deNovoAccFlag': 1,
                                                                               'refZScore': -9.017236678144027})
    @mock.patch('calcVarPriors.getPriorProbDeNovoDonorSNS', return_value = {'refMaxEntScanScore': -14.15,
                                                                            'altMaxEntScanScore': -5.87,
                                                                            'enigmaClass': 'class_2',
                                                                            'altZScore': -5.928774792466758,
                                                                            'priorProb': 0.02,
                                                                            'deNovoDonorFlag': 1,
                                                                            'refZScore': -9.483955273593583})
    @mock.patch('calcVarPriors.getPriorProbProteinSNS', return_value = {'enigmaClass': 'class_2',
                                                                        'priorProb': 0.02})
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "splice_region_variant")
    def test_getPriorProbSpliceAcceptorSNSNoDeNovoAccBRCA2(self, varInSpliceRegion, getVarType, varInExon,
                                                           getPriorProbRefSpliceAccceptorSNS, getPriorProbDeNovoAcceptorSNS,
                                                           getPriorProbDeNovoDonorSNS, getPriorProbProteinSNS, getVarConsequences):
        '''Tests that applicable prior for a variant in a reference splice site is assigned correctly (no de novo splicing)'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Chr"] = "13"
        self.variant["HGVS_cDNA"] = "c.7008C>G"
        self.variant["Pos"] = self.variant["Hg38_Start"] = self.variant["Hg38_End"] = "32354861"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "G"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries, variantData)
        # checks that variant splice site flag is assigned correctly
        self.assertEquals(priorProb["spliceSite"], 1)
        # checks that variant is flagged as a de novo splice donor and acceptor
        self.assertEquals(priorProb["deNovoDonorFlag"], 1)
        self.assertEquals(priorProb["deNovoAccFlag"], 1)
        # checks that prior prob and enigma class are appropriate based on applicable prior
        self.assertEquals(priorProb["applicablePrior"], priorProbs["low"])
        self.assertEquals(priorProb["applicableEnigmaClass"], enigmaClasses["class2"])
        # checks that protein prior prob, ref prior prob, and de novo prior probs are set appropriately
        self.assertEquals(priorProb["proteinPrior"], priorProbs["deNovoLow"])
        self.assertEquals(priorProb["refAccPrior"], priorProbs["low"])
        self.assertEquals(priorProb["deNovoAccPrior"], "N/A")
        self.assertEquals(priorProb["deNovoDonorPrior"], priorProbs["deNovoLow"])
        # checks that a score is present for reference acceptor value and de novo donor/acceptor value
        self.assertNotEquals(priorProb["refRefAccZ"], "-")
        self.assertNotEquals(priorProb["refDeNovoAccZ"], "-")
        self.assertNotEquals(priorProb["refDeNovoDonorZ"], "-")
        # checks that scores are not present for ref splice donor site or de novo splice acceptor sites
        self.assertEquals(priorProb["refRefDonorZ"], "-")
        # checks that splice rescue, splice flag, and frameshift are all equal to zero
        self.assertEquals(priorProb["spliceRescue"], 0)
        self.assertEquals(priorProb["spliceFlag"], 0)
        self.assertEquals(priorProb["frameshift"], 0)
        
    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = False)
    @mock.patch('calcVarPriors.getPriorProbRefSpliceAcceptorSNS', return_value = {'refMaxEntScanScore': 4.57,
                                                                                  'altMaxEntScanScore': 2.42,
                                                                                  'enigmaClass': 'class_4',
                                                                                  'altZScore': -2.286639792272834,
                                                                                  'priorProb': 0.97,
                                                                                  'spliceSite': 1,
                                                                                  'refZScore': -1.4031975880833913})
    @mock.patch('calcVarPriors.getPriorProbDeNovoAcceptorSNS', return_value = {'refMaxEntScanScore': -1.89,
                                                                               'altMaxEntScanScore': 6.71,
                                                                               'enigmaClass': 'N/A',
                                                                               'altZScore': -0.5238644174018071,
                                                                               'priorProb': 'N/A',
                                                                               'deNovoAccFlag': 1,
                                                                               'refZScore': -4.057633234159576})
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "intron_variant")
    def test_getPriorProbSpliceAcceptorSNSWithDeNovoBRCA1(self, varInSpliceRegion, getVarType, varInExon,
                                                          getPriorProbRefSpliceAccceptorSNS, getPriorProbDeNovoAcceptorSNS,
                                                          getVarConsequences):
        '''Tests that applicable prior for a variant in a reference splice site is assigned correctly (with de novo splicing)'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Chr"] = 17
        self.variant["HGVS_cDNA"] = "c.4186-16t>G"
        self.variant["Pos"] = self.variant["Hg38_Start"] = self.variant["Hg38_End"] = "43082591"
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "C"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries, variantData)
        # checks that variant splice site flag is assigned correctly
        self.assertEquals(priorProb["spliceSite"], 1)
        # checks that variant is flagged as a de novo splice acceptor
        self.assertEquals(priorProb["deNovoAccFlag"], 1)
        # checks that prior prob and enigma class are appropriate based on applicable prior
        self.assertEquals(priorProb["applicablePrior"], priorProbs["high"])
        self.assertEquals(priorProb["applicableEnigmaClass"], enigmaClasses["class4"])
        # checks that protein prior prob, ref prior prob, and de novo prior prob are set appropriately
        self.assertEquals(priorProb["proteinPrior"], priorProbs["NA"])
        self.assertEquals(priorProb["refAccPrior"], priorProbs["high"])
        self.assertEquals(priorProb["deNovoAccPrior"], priorProbs["NA"])
        # checks that a score is present for reference acceptor score and de novo splice acceptor score
        self.assertNotEquals(priorProb["altDeNovoAccMES"], "-")
        self.assertNotEquals(priorProb["refRefAccZ"], "-")
        # checks that scores are not present for ref splice donor site or de novo splice donor sites
        self.assertEquals(priorProb["altDeNovoDonorMES"], "-")
        self.assertEquals(priorProb["refRefDonorZ"], "-")
        # checks that splice rescue, splice flag, and frameshift are all equal to zero
        self.assertEquals(priorProb["spliceRescue"], 0)
        self.assertEquals(priorProb["spliceFlag"], 0)
        self.assertEquals(priorProb["frameshift"], 0)

    @mock.patch('calcVarPriors.varInSpliceRegion', return_value = True)
    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    @mock.patch('calcVarPriors.varInExon', return_value = True)
    @mock.patch('calcVarPriors.getPriorProbRefSpliceAcceptorSNS', return_value = {'refMaxEntScanScore': 3.71,
                                                                                  'altMaxEntScanScore': 3.02,
                                                                                  'enigmaClass': 'class_3',
                                                                                  'altZScore': -2.0400977818013613,
                                                                                  'priorProb': 0.34,
                                                                                  'spliceSite': 1,
                                                                                  'refZScore': -1.7565744697591685})
    @mock.patch('calcVarPriors.getPriorProbDeNovoDonorSNS', return_value = {'refMaxEntScanScore': 0.48,
                                                                            'altMaxEntScanScore': -1.52,
                                                                            'enigmaClass': 'N/A',
                                                                            'altZScore': -4.061016931005201,
                                                                            'priorProb': 'N/A',
                                                                            'deNovoDonorFlag': 0,
                                                                            'refZScore': -3.2022776843562095})
    @mock.patch('calcVarPriors.getPriorProbDeNovoAcceptorSNS', return_value = {'refMaxEntScanScore': 4.78,
                                                                               'altMaxEntScanScore': 3.43,
                                                                               'enigmaClass': 'N/A',
                                                                               'altZScore': -1.8716274079791886,
                                                                               'priorProb': 'N/A',
                                                                               'deNovoAccFlag': 0,
                                                                               'refZScore': -1.3169078844183761})
    @mock.patch('calcVarPriors.getPriorProbProteinSNS', return_value = {'enigmaClass': 'class_5',
                                                                        'priorProb': 0.99})
    @mock.patch('calcVarPriors.getVarConsequences', return_value = "stop_gained")
    @mock.patch('calcVarPriors.getPriorProbSpliceRescueNonsenseSNS', return_value = {'spliceFlag': 0,
                                                                                     'enigmaClass': 'class_4',
                                                                                     'frameshift': 0,
                                                                                     'priorProb': 0.97,
                                                                                     'spliceRescue': 0})
    def test_getPriorProbSpliceAcceptorSNSNonsenseBRCA1(self, varInSpliceRegion, getVarType, varInExon,
                                                        getPriorProbRefSpliceAcceptorSNS, getPriorProbDeNovoDonorSNS,
                                                        getPriorProbDeNovoAcceptorSNS, getPriorProbProteinSNS,
                                                        getVarConsequences, getPriorProbSpliceRescueNonsenseSNS):
        '''Tests function with nonsense variant in reference acceptor site'''
        boundaries = "enigma"
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Chr"] = "17"
        self.variant["HGVS_cDNA"] = "c.442C>T"
        self.variant["Pos"] = self.variant["Hg38_Start"] = self.variant["Hg38_End"] = "43099880"
        self.variant["Ref"] = "G"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbSpliceAcceptorSNS(self.variant, boundaries, variantData)
        # checks that variant splice site flag is assigned correctly
        self.assertEquals(priorProb["spliceSite"], 1)
        # checks that variant is not flagged as a de novo splice donor or acceptor
        self.assertEquals(priorProb["deNovoDonorFlag"], 0)
        self.assertEquals(priorProb["deNovoAccFlag"], 0)
        # checks that prior prob and enigma class are appropriate based on applicable prior
        self.assertEquals(priorProb["applicablePrior"], priorProbs["high"])
        self.assertEquals(priorProb["applicableEnigmaClass"], enigmaClasses["class4"])
        # checks that protein prior prob, ref prior prob, and de novo prior prob are set appropriately
        self.assertEquals(priorProb["proteinPrior"], priorProbs["pathogenic"])
        self.assertEquals(priorProb["refAccPrior"], priorProbs["moderate"])
        self.assertEquals(priorProb["deNovoDonorPrior"], priorProbs["NA"])
        self.assertEquals(priorProb["deNovoAccPrior"], priorProbs["NA"])
        # checks that a score is NOT present for reference donor score and is present for de novo splice donor score
        self.assertNotEquals(priorProb["altDeNovoDonorZ"], "-")
        self.assertEquals(priorProb["refRefDonorMES"], "-")
        # checks that scores are present for ref splice acceptor site or de novo splice acceptor sites
        self.assertNotEquals(priorProb["altDeNovoAccZ"], "-")
        self.assertNotEquals(priorProb["refRefAccMES"], "-")
        # checks that splice rescue, splice flag, and frameshift are equal to zero
        self.assertEquals(priorProb["spliceRescue"], 0)
        self.assertEquals(priorProb["spliceFlag"], 0)
        self.assertEquals(priorProb["frameshift"], 0)

    @mock.patch('calcVarPriors.getVarType', return_value = varTypes["sub"])
    def test_getPriorProbProteinSNS(self, getVarType):
        '''Tests that function parses data from variantData correctly and returns correct prior prob/class'''
        # checks for BRCA1 variant
        self.variant["Gene_Symbol"] = "BRCA1"
        self.variant["HGVS_cDNA"] = "c.592A>T"
        self.variant["Ref"] = "T"
        self.variant["Alt"] = "A"
        priorProb = calcVarPriors.getPriorProbProteinSNS(self.variant, variantData)
        self.assertEquals(priorProb["priorProb"], priorProbs["proteinMod"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

        # checks for BRCA2 variant
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["HGVS_cDNA"] = "c.620C>T"
        self.variant["Ref"] = "C"
        self.variant["Alt"] = "T"
        priorProb = calcVarPriors.getPriorProbProteinSNS(self.variant, variantData)
        self.assertEquals(priorProb["priorProb"], priorProbs["proteinHigh"])
        self.assertEquals(priorProb["enigmaClass"], enigmaClasses["class3"])

    @mock.patch('calcMaxEntScanMeanStd.fetch_gene_coordinates', return_value = transcriptDataBRCA2)
    def test_getVarDict(self, fetch_gene_coordinates):
        '''
        Tests that: 
        1. Variant information is being parsed correctly
        '''
        boundaries = "enigma"
        varDict = calcVarPriors.getVarDict(self.variant, boundaries)
        self.assertEquals(varDict["varHGVScDNA"], self.variant["HGVS_cDNA"])
        self.assertEquals(varDict["varChrom"], self.variant["Chr"])
        self.assertEquals(varDict["varGene"], self.variant["Gene_Symbol"])
        self.assertEquals(varDict["varGenCoordinate"], self.variant["Pos"])
