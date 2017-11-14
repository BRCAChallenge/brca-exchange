import unittest
import calcVarPriors

class test_calcVarPriors(unittest.TestCase):

    def setUp(self):

        self.variant = {"Chr":"13",
                        "Pos":"32314943",
                        "Ref":"A",
                        "Alt":"G",
                        "Gene_Symbol":"BRCA2",
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
