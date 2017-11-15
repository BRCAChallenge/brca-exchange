import unittest
import calcVarPriors

class test_calcVarPriors(unittest.TestCase):

    def setUp(self):

        self.variant = {"Chr":"13",
                        "Pos":"32314943",
                        "Ref":"A",
                        "Alt":"G",
                        "Gene_Symbol":"BRCA2",
                        "Reference_Sequence": "NM_000059.3",
                        "pyhgvs_cDNA":"NM_000059.3:c.-764A>G"}
                           
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

        
    def test_getVarType(self):
        '''
        Tests that variant type is set correctly to substitution, deletion, insertion, or delins based on variant "Ref" and "Alt" values
        '''
        self.variant["Ref"] = "A"
        self.variant["Alt"] = "T"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, "substitution")

        self.variant["Ref"] = "A"
        self.variant["Alt"] = "AAA"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, "insertion")

        self.variant["Ref"] = "AGT"
        self.variant["Alt"] = "A"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, "deletion")

        self.variant["Ref"] = "AG"
        self.variant["Alt"] = "AGTA"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, "delins")

        self.variant["Ref"] = "AGTA"
        self.variant["Alt"] = "AG"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, "delins")

        self.variant["Ref"] = "AG"
        self.variant["Alt"] = "GT"
        varType = calcVarPriors.getVarType(self.variant)
        self.assertEquals(varType, "delins")

    def test_getExonBoundaries(self):
        '''
        Tests that getExonBoundaries returns the correct exon boundaries for a gene on the positive strand (BRCA2) and a gene on the negative strand (BRCA1) 
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"    # BRCA1
        varExons = calcVarPriors.getExonBoundaries(self.variant["Reference_Sequence"])
        # check correct number of exons
        numExons = 23
        self.assertEquals(len(varExons), numExons)
        # because '-' strand gene, checks that exonEnd > exonStart
        for exon in varExons.keys():
            exonBounds = varExons[exon]
            exonStart = exonBounds['exonStart']    
            exonEnd = exonBounds['exonEnd'] 
            self.assertGreater(exonStart, exonEnd)
            exonNum = int(exon[4:])      # gets current exon number
            nextExonNum = "exon" + str(exonNum + 1)  # gets key for next exon number
            if nextExonNum <= numExons:
                nextExonStart = varExons[nextExonNum]['exonStart']
                # checks to make sure that next exon does not start until after current exon ends
                self.assertGreater(exonEnd, nextExonStart)

        self.variant["Reference_Sequence"] = "NM_000059.3"    # BRCA2
        varExons = calcVarPriors.getExonBoundaries(self.variant["Reference_Sequence"])
        # check correct number of exons
        numExons = 27
        self.assertEquals(len(varExons), numExons)
        # because '+' strand gene, checks that exonEnd > exonStart
        for exon in varExons.keys():
            exonBounds = varExons[exon]
            exonStart = exonBounds['exonStart']
            exonEnd = exonBounds['exonEnd']
            self.assertGreater(exonEnd, exonStart)
            exonNum = int(exon[4:])      # get current exon number
            nextExonNum = "exon" + str(exonNum + 1)   # gets key for next exon number
            if nextExonNum <= numExons:
                nextExonStart = varExons[nextExonNum]['exonStart']
                # checks to make sure that next exon does not start until after current exon ends
                self.assertGreater(nextExonStart, exonEnd)

    def test_varOutsideBoundaries(self):
        '''Tests that varOutsideBoundaries correctly classifies variants based on genomic posiiton'''
        # checks BRCA1 variant outside transcript boundaries
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43044274"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertTrue(varOutBounds)

        # checks BRCA1 variant inside transcript boundaries
        self.variant["Pos"] = "43070957"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertFalse(varOutBounds)

        # checks BRCA2 variant outside transcript boundaries
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32315477"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertTrue(varOutBounds)

        #checks BRCA2 variant inside transcript boundaries
        self.variant["Pos"] = "32326500"
        varOutBounds = calcVarPriors.varOutsideBoundaries(self.variant)
        self.assertFalse(varOutBounds)

    def test_varInExon(self):    
        '''Tests that variants are correctly classified as in an exon or outside of an exon'''
        # checks BRCA1 variant in an exon
        self.variant["Reference_Sequence"] = "NM_007294.3"
        self.variant["Pos"] = "43097273"
        varExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(varExon)

        # checks BRCA1 variant outside an exon
        self.variant["Pos"] = "43047733"
        varExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(varExon)
        
        # checks BRCA2 variant in an exon
        self.variant["Reference_Sequence"] = "NM_000059.3"
        self.variant["Pos"] = "32354890"
        varExon = calcVarPriors.varInExon(self.variant)
        self.assertTrue(varExon)
        
        # checks BRCA2 variant outside an exon
        self.variant["Pos"] = "32379555"
        varExon = calcVarPriors.varInExon(self.variant)
        self.assertFalse(varExon)
                
    def test_getVarLocation(self):
        '''
        Tests that:
        1. Variant location is set correctly for genomic position outside transcript boundaries
        2. Variant location is set correctly for genomic position in exon
        3. Variant location is set correclty for genomic position in intron
        '''
        self.variant["Reference_Sequence"] = "NM_007294.3"
        # position before txn start site for BRCA1
        self.variant["Pos"] = "43044274"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        self.assertEquals(varLoc, "outsideBoundaries")

        # position in exon for BRCA1
        self.variant["Pos"] = "43070957"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        self.assertEquals(varLoc, "inExon")
        
        # position in intron for BRCA1
        self.variant["Pos"] = "43091062"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        self.assertEquals(varLoc, "inIntron")

        self.variant["Reference_Sequence"] = "NM_000059.3"
        # position before txn start site for BRCA2
        self.variant["Pos"] = "32315477"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        self.assertEquals(varLoc, "outsideBoundaries")

        # position in exon for BRCA2
        self.variant["Pos"] = "32326500"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        self.assertEquals(varLoc, "inExon")
        
        # position in intron for BRCA2
        self.variant["Pos"] = "32330915"
        varLoc = calcVarPriors.getVarLocation(self.variant)
        self.assertEquals(varLoc, "inIntron")
        

    def test_getVarDict(self):
        '''
        Tests that: 
        1. Variant information is being parsed correctly
        2. Variant strand is set correctly based on variant gene
        '''

        varDict = calcVarPriors.getVarDict(self.variant)
        self.assertEquals(varDict["varHGVScDNA"], self.variant["pyhgvs_cDNA"])
        self.assertEquals(varDict["varChrom"], self.variant["Chr"])
        self.assertEquals(varDict["varGene"], self.variant["Gene_Symbol"])
        self.assertEquals(varDict["varGenCoordinate"], self.variant["Pos"])
        
        self.variant["Gene_Symbol"] = "BRCA1"
        varDict = calcVarPriors.getVarDict(self.variant)
        self.assertEquals(varDict["varStrand"], "-")

        self.variant["Gene_Symbol"] = "BRCA2"
        varDict = calcVarPriors.getVarDict(self.variant)
        self.assertEquals(varDict["varStrand"], "+")
