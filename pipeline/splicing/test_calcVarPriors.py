import unittest
import calcVarPriors

class test_calcVarPriors(unittest.TestCase):

    def setUp(self):

        self.variant = {"Chr": "13",
                        "Pos": "32314943",
                        "Ref": "A",
                        "Alt": "G",
                        "Gene_Symbol": "BRCA2",
                        "Reference_Sequence": "NM_000059.3",
                        "pyhgvs_cDNA": "NM_000059.3:c.-764A>G"}

        self.strand = {"minus": "-",
                       "plus": "+"}

        self.varTypes = {"sub": "substitution",
                         "ins": "insertion",
                         "del": "deletion",
                         "delins": "delins"}

        self.numExons = {"BRCA1": 23,
                         "BRCA2": 27}

        self.exonDonorBoundsBRCA1 = {"exon16": {"donorStart": 43070930,
                                                "donorEnd": 43070922}}

        self.exonDonorBoundsBRCA2 = {"exon15": {"donorStart": 32356607,
                                                "donorEnd": 32356615}}

        self.exonAcceptorBoundsBRCA1 = {"exon21": {"acceptorStart": 43051137,
                                                   "acceptorEnd": 43051115}}

        self.exonAcceptorBoundsBRCA2 = {"exon20": {"acceptorStart": 32370936,
                                                   "acceptorEnd": 32370958}}
                           
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
        self.assertEquals(varStrand, self.strand["minus"])

        self.variant["Gene_Symbol"] = "BRCA2"
        varStrand = calcVarPriors.getVarStrand(self.variant)
        self.assertEquals(varStrand, self.strand["plus"])

        
    def test_getVarType(self):
        '''
        Tests that variant type is set correctly to substitution, deletion, insertion, or delins based on variant "Ref" and "Alt" values
        '''
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, self.varTypes["sub"])

        self.variant["Ref"] = "A"
        self.variant["Alt"] = "AAA"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, self.varTypes["ins"])

        self.variant["Ref"] = "AGT"
        self.variant["Alt"] = "A"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, self.varTypes["del"])

        self.variant["Ref"] = "AG"
        self.variant["Alt"] = "AGTA"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, self.varTypes["delins"])

        self.variant["Ref"] = "AGTA"
        self.variant["Alt"] = "AG"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, self.varTypes["delins"])

        self.variant["Ref"] = "AG"
        self.variant["Alt"] = "GT"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, self.varTypes["delins"])

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


    def test_varOutsideBoundaries(self):
        '''Tests that variant outside/inside transcript boundaries are correctly identified'''
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

    def test_getExonBoundaries(self):
        '''
        Tests that:
        1. Exon boundaries are set correctly for:
        gene on minus strand (BRCA1) and gene on plus strand (BRCA2)
            - checks that exon does not start before it ends
            - next exon does not start before current exon ends
        2. length of varExons matches number of exons for a specific gene
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        varExons = calcVarPriors.getExonBoundaries(self.variant)
        # check correct number of exons
        self.assertEquals(len(varExons), self.numExons["BRCA1"])
        # because '-' strand gene, checks that exonStart > exonEnd
        for exon in varExons.keys():
            exonBounds = varExons[exon]
            self.assertGreater(exonBounds["exonStart"], exonBounds["exonEnd"])
            currentExonNum = int(exon[4:])
            nextExonNum = str(currentExonNum + 1)
            nextExonKey = "exon" + nextExonNum
            if nextExonNum <= self.numExons["BRCA1"]:
                nextExonStart = varExons[nextExonKey]["exonStart"]
                # checks that next exon does not start before current exon ends
                self.assertGreater(exonBounds["exonEnd"], nextExonStart)

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        varExons = calcVarPriors.getExonBoundaries(self.variant)
        # check correct number of exons
        self.assertEquals(len(varExons), self.numExons["BRCA2"])
        # because '+' strand gene, checks that exonEnd > exonStart
        for exon in varExons.keys():
            exonBounds = varExons[exon]
            self.assertGreater(exonBounds["exonEnd"], exonBounds["exonStart"])
            currentExonNum = int(exon[4:])
            nextExonNum = str(currentExonNum + 1)
            nextExonKey = "exon" + nextExonNum
            if nextExonNum <= self.numExons["BRCA2"]:
                nextExonStart = varExons[nextExonKey]["exonStart"]
                # cehcks that next exon does not start before current exon ends
                self.assertGreater(nextExonStart, exonBounds["exonEnd"])


    def test_getRefSpliceDonorBoundaries(self):
        '''
        Tests that splice donor boundaries are set correctly for reference transcript and strand
        Uses example exon boundaries set in setUp function
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        spliceDonorBounds = calcVarPriors.getRefSpliceDonorBoundaries(self.variant)
        # to find exon specified in setUp function
        exon = self.exonDonorBoundsBRCA1.keys()[0]
        self.assertEquals(self.exonDonorBoundsBRCA1[exon]["donorStart"],
                          spliceDonorBounds[exon]["donorStart"])
        self.assertEquals(self.exonDonorBoundsBRCA1[exon]["donorEnd"],
                          spliceDonorBounds[exon]["donorEnd"])

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        spliceDonorBounds = calcVarPriors.getRefSpliceDonorBoundaries(self.variant)
        # to find exon specified in setUp function
        exon = self.exonDonorBoundsBRCA2.keys()[0]
        self.assertEquals(self.exonDonorBoundsBRCA2[exon]["donorStart"],
                          spliceDonorBounds[exon]["donorStart"])
        self.assertEquals(self.exonDonorBoundsBRCA2[exon]["donorEnd"],
                          spliceDonorBounds[exon]["donorEnd"])


    def test_getRefSpliceAcceptorBoundaries(self):
        '''
        Tests that splice acceptor boundaries are set correctly for refernce transcript and strand
        Uses example exon boundaries in setUp function
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        spliceAcceptorBounds = calcVarPriors.getRefSpliceAcceptorBoundaries(self.variant)
        # to find exon specified in setUp function
        exon = self.exonAcceptorBoundsBRCA1.keys()[0]
        self.assertEquals(self.exonAcceptorBoundsBRCA1[exon]["acceptorStart"],
                          spliceAcceptorBounds[exon]["acceptorStart"])
        self.assertEquals(self.exonAcceptorBoundsBRCA1[exon]["acceptorEnd"],
                          spliceAcceptorBounds[exon]["acceptorEnd"])

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        spliceAcceptorBounds = calcVarPriors.getRefSpliceAcceptorBoundaries(self.variant)
        # to find exon specified in setUp function
        exon = self.exonAcceptorBoundsBRCA2.keys()[0]
        self.assertEquals(self.exonAcceptorBoundsBRCA2[exon]["acceptorStart"],
                          spliceAcceptorBounds[exon]["acceptorStart"])
        self.assertEquals(self.exonAcceptorBoundsBRCA2[exon]["acceptorEnd"],
                          spliceAcceptorBounds[exon]["acceptorEnd"])


    def test_varInExon(self):
        '''Tests that variant is correctly identified as inside or outside an exon'''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"

        # checks BRCA1 variant inside an exon
        self.variant["Pos"] = "43049176"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # checks BRCA1 variant outside an exon
        self.variant["Pos"] = "43045827"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"

        # checks BRCA2 variant inside an exon
        self.variant["Pos"] = "32398201"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(inExon)

        # cehcks BRCA2 variant outside an exon
        self.variant["Pos"] = "32396873"
        inExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(inExon)

    def test_varInSpliceDonor(self):
        '''
        Tests that variant is correctly identified as in or NOT in a splice donor region 
        for multiple positions across multiple exons
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"

        #checks that 7th base in intron is NOT counted as in splice donor
        self.variant["Pos"] = "43063326"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertFalse(inSpliceDonor)
        
        # checks that first base in intron counted as in splice donor
        self.variant["Pos"] =  "43097243"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in intron counted as in splice donor
        self.variant["Pos"] = "43095842"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that last base in intron counted as in splice donor
        self.variant["Pos"] = "43091429"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that first base in exon counted as in splice donor
        self.variant["Pos"] = "43090946"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in exon counted as in splice donor
        self.variant["Pos"] = "43082405"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that last base in exon counted as in splice donor
        self.variant["Pos"] = "43076488"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that 4th to last base in exon is NOT counted as in splice donor
        self.variant["Pos"] = "43057055"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertFalse(inSpliceDonor)

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        self.variant["Chr"] = "13"

        # checks that 7th base in intron is NOT counted as in splice donor
        self.variant["Pos"] = "32363540"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertFalse(inSpliceDonor)

        # checks that first base in intron counted as in splice donor
        self.variant["Pos"] = "32329493"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in intron counted as in splice donor
        self.variant["Pos"] = "32331033"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that last base in intron counted as in splice donor
        self.variant["Pos"] = "32333393"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that first base in exon counted as in splice donor
        self.variant["Pos"] = "32341194"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that middle base in exon counted as in splice donor
        self.variant["Pos"] = "32344652"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that last base in exon counted as in splice donor
        self.variant["Pos"] = "32346896"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertTrue(inSpliceDonor)

        # checks that 4th to last base in exon is NOT counted as in splice donor
        self.variant["Pos"] = "32370554"
        inSpliceDonor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertFalse(inSpliceDonor)


    def test_varInSpliceAcceptor(self):
        '''
        Tests that variant is correctly identified as in or NOT in a splice acceptor region
        uses multiple positions across multiple exons
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"

        # checks that -21st base in intron is NOT counted as in splice acceptor
        self.variant["Pos"] = "43067716"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertFalse(inSpliceAcceptor)
        
        # checks that first base in intron counted as in splice acceptor
        self.variant["Pos"] = "43124135"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in intron counted as in splice acceptor
        self.variant["Pos"] = "43115787"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in intron counted as in splice acceptor
        self.variant["Pos"] = "43106534"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that first base in exon counted as in splice acceptor
        self.variant["Pos"] = "43104956"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in exon counted as in splice acceptor
        self.variant["Pos"] = "43104260"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in exon counted as in splice acceptor
        self.variant["Pos"] = "43099878"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that 4th base in exon is NOT counted as in splice acceptor
        self.variant["Pos"] = "43063948"
        inSpliceAcceptor = calcVarPriors.varInSpliceDonor(self.variant)
        self.assertFalse(inSpliceAcceptor)

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"

        # checks that -21st base in intron is NOT counted as in splice acceptor
        self.variant["Pos"] = "32357721"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertFalse(inSpliceAcceptor)

        # checks that first base in intron counted as in splice acceptor
        self.variant["Pos"] = "32316402"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in intron counted as in splice acceptor
        self.variant["Pos"] = "32319069"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in intron counted as in splice acceptor
        self.variant["Pos"] = "32325075"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that first base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326101"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that middle base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326243"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that last base in exon counted as in splice acceptor
        self.variant["Pos"] = "32326501"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertTrue(inSpliceAcceptor)

        # checks that 4th base in exon is NOT counted as in splice acceptor
        self.variant["Pos"] = "32362526"
        inSpliceAcceptor = calcVarPriors.varInSpliceAcceptor(self.variant)
        self.assertFalse(inSpliceAcceptor)
        
                
    def test_getVarLocation(self):
        '''
        Tests that:
        1. Variant location is set correctly for genomic position outside transcript boundaries
        2. Variant location is set correctly for genomic position in exon
        3. Variant location is set correclty for genomic position in intron
        '''
        # TO DO - implement tests for varLocation once varLocation changed to use Ensembl API
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Gene_Symbol"] = "BRCA1"
        # position before txn start site for BRCA1
        self.variant["Pos"] = "43044274"
        varLoc = calcVarPriors.getVarLocation(self.variant)

        # position in exon for BRCA1
        self.variant["Pos"] = "43070957"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        
        # position in intron for BRCA1
        self.variant["Pos"] = "43106443"
        varLoc = calcVarPriors.getVarLocation(self.variant)

        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Gene_Symbol"] = "BRCA2"
        # position before txn start site for BRCA2
        self.variant["Pos"] = "32315477"
        varLoc = calcVarPriors.getVarLocation(self.variant)

        # position in exon for BRCA2
        self.variant["Pos"] = "32326500"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        
        # position in intron for BRCA2
        self.variant["Pos"] = "32357952"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        

    def test_getVarDict(self):
        '''
        Tests that: 
        1. Variant information is being parsed correctly
        '''

        varDict = calcVarPriors.getVarDict(self.variant)
        self.assertEquals(varDict["varHGVScDNA"], self.variant["pyhgvs_cDNA"])
        self.assertEquals(varDict["varChrom"], self.variant["Chr"])
        self.assertEquals(varDict["varGene"], self.variant["Gene_Symbol"])
        self.assertEquals(varDict["varGenCoordinate"], self.variant["Pos"])
        
