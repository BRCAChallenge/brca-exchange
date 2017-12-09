import pytest
import unittest
import checkVarHGVS

class testCheckVarHGVS(unittest.TestCase):

    def setUp(self):

        self.varBRCA1 = {"Ref": "C",
                         "Alt": "T",
                         "pyhgvs_cDNA": "NM_007294.3:c.5485G>A",
                         "Gene_Symbol": "BRCA1"}
        self.varBRCA2 = {"Ref": "A",
                         "Alt": "T",
                         "pyhgvs_cDNA": "NM_000059.3:c.66A>T",
                         "Gene_Symbol": "BRCA2"}

        self.sequence = "AGT"

    def test_getRevComp(self):
        '''Test funciton getRevComp'''
        # checks that getRevComp works for single nucleotide
        self.sequence = "A"
        revComp = checkVarHGVS.getRevComp(self.sequence)
        self.assertEquals(revComp, "T")

        # checks that getRevComp works for multiple nucleotides
        self.sequence = "AGT"
        revComp = checkVarHGVS.getRevComp(self.sequence)
        self.assertEquals(revComp, "ACT")

        # checks that getRevComp works for multiple nucleotides with unspectiifed nucleotides
        self.sequence = "AGTNT"
        revComp = checkVarHGVS.getRevComp(self.sequence)
        self.assertEquals(revComp, "ANACT") 

    def test_getVarInfo(self):
        '''Tests function getVarInfo'''
        # checks that variant is being parsed correctly
        varInfo = checkVarHGVS.getVarInfo(self.varBRCA1)
        self.assertEquals(varInfo["varcDNAHGVS"], self.varBRCA1["pyhgvs_cDNA"])
        self.assertEquals(varInfo["varGenRef"], self.varBRCA1["Ref"])
        self.assertEquals(varInfo["varGenAlt"], self.varBRCA1["Alt"])
        self.assertEquals(varInfo["varGene"], self.varBRCA1["Gene_Symbol"])

        # checks that cDNARef and cDNAAlt are set correctly for BRCA1 (reverse complement)
        self.varBRCA1["Gene_Symbol"] = "BRCA1"
        varInfo = checkVarHGVS.getVarInfo(self.varBRCA1)
        self.assertEquals(varInfo["cDNARef"], self.varBRCA1["Ref"])
        self.assertEquals(varInfo["cDNAAlt"], self.varBRCA1["Alt"])

        # checks that cDNARef and cDNAAlt are set correctly for BRCA2
        self.varBRCA2["Gene_Symbol"] = "BRCA2"
        varInfo = checkVarHGVS.getVarInfo(self.varBRCA2)
        self.assertEquals(varInfo["cDNARef"], self.varBRCA2["Ref"])
        self.assertEquals(varInfo["cDNAAlt"], self.varBRCA2["Alt"])

    def test_checkVarHGVSSub(self):
        '''Tests checkVarHGVS for single nucleotide substitutions'''
        # checks that checkVarHGVS works for substitution in BRCA1
        self.varBRCA1["Ref"] = "C"
        self.varBRCA1["Alt"] = "T"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.5485G>A"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "T"
        self.varBRCA1["Alt"] = "G"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)
        
        # checks that checkVarHGVS works for substitution in BRCA2
        self.varBRCA2["Ref"] = "A"
        self.varBRCA2["Alt"] = "T"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.66A>T"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "T"
        self.varBRCA2["Alt"] = "G"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)

    def test_checkVarHGVSDelSingle(self):
        '''Tests checkVarHGVS for single nucleotide deletion variants'''
        # checks that checkVarHGVS works for single nucleotide deletion in BRCA1
        self.varBRCA1["Ref"] = "GA"
        self.varBRCA1["Alt"] = "G"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.*179delT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "TG"
        self.varBRCA1["Alt"] = "T"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)

        # checks that checkVarHGVS works for single nucleotide deletion in BRCA2
        self.varBRCA2["Ref"] = "TT"
        self.varBRCA2["Alt"] = "T"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.-39-5delT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "GG"
        self.varBRCA2["Alt"] = "G"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)

    def test_checkVarHGVSDelMult(self):
        '''Tests checkVarHGVS for multiple nucleotide deletions'''
        # checks that checkVarHGVS works for a multiple nucleotide deletion in BRCA1
        self.varBRCA1["Ref"] = "ATAAAG"
        self.varBRCA1["Alt"] = "A"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.-19-85_-19-81delCTTTA"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "TATTA"
        self.varBRCA1["Alt"] = "T"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)
        
        # checks that checkVarHGVS works for a multiple nucleotide deletion in BRCA2
        self.varBRCA2["Ref"] = "GCG"
        self.varBRCA2["Alt"] = "G"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.-42_-41delCG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "TAGA"
        self.varBRCA2["Alt"] = "T"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)

    def test_checkVarHGVSInsSingle(self):
        '''Tests checkVarHGVS for single nucleotide insertions'''
        # checks that checkVarHGVS works for a single nucleotide insertion in BRCA1
        self.varBRCA1["Ref"] = "C"
        self.varBRCA1["Alt"] = "CG"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.1017_1018insC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "T"
        self.varBRCA1["Alt"] = "TT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)

        # checks that checkVarHGVS works for a single nucleotide insertion in BRCA2
        self.varBRCA2["Ref"] = "C"
        self.varBRCA2["Alt"] = "CA"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.1584_1585insA"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "G"
        self.varBRCA2["Alt"] = "GG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)

    def test_checkVarHGVSInsMult(self):
        '''Tests checkVarHGVS for multiple nucleotide insertions'''
        # checks that checkVarHGVS works for a multiple nucleotide insertion in BRCA1
        self.varBRCA1["Ref"] = "A"
        self.varBRCA1["Alt"] = "ACTTTC"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.1685_1686insGAAAG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "T"
        self.varBRCA1["Alt"] = "TATTAC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)
        
        # checks that checkVarHGVS works for a multiple nucleotide insertion in BRCA2
        self.varBRCA2["Ref"] = "A"
        self.varBRCA2["Alt"] = "ATTAG"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.1190_1191insTTAG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "A"
        self.varBRCA2["Alt"] = "ATTAC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)
        
    def test_checkVarHGVSDupSingle(self):
        '''Tests checkVarHGVS for single nucleotide duplications'''
        # checks that checkVarHGVS works for a single nucleotide duplication in BRCA1
        self.varBRCA1["Ref"] = "T"
        self.varBRCA1["Alt"] = "TT"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.1008dupA"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "G"
        self.varBRCA1["Alt"] = "GG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)

        # checks that checkVarHGVS works for a single nucleotide duplication in BRCA2

        self.varBRCA2["Ref"] = "A"
        self.varBRCA2["Alt"] = "AC"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.103dupC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "T"
        self.varBRCA2["Alt"] = "TG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)

        
    def test_checkVarHGVSDupMult(self):
        '''Tests checkVarHGVS for multiple nucleotide duplications'''
        # checks that checkVarHGVS works for a multiple nucleotide duplication in BRCA1
        self.varBRCA1["Ref"] = "A"
        self.varBRCA1["Alt"] = "ATA"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.-19-22_-19-21dupAT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "G"
        self.varBRCA1["Alt"] = "GCG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)
        
        # checks that checkVarHGVS works for a multiple nucleotide duplication in BRCA2
        self.varBRCA2["Ref"] = "G"
        self.varBRCA2["Alt"] = "GTACC"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.1163_1166dupTACC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "T"
        self.varBRCA2["Alt"] = "TGGAC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)


    def test_checkVarHGVSDelinsEqual(self):
        '''Tests checkVarHGVS for delins variants with equal number of nucleotides deleted and inserted'''
        
        # checks that checkVarHGVS works for an equal nucleotide delins in BRCA1
        self.varBRCA1["Ref"] = "GA"
        self.varBRCA1["Alt"] = "AG"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.441+51_441+52delTCinsCT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "GG"
        self.varBRCA1["Alt"] = "TT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)
        
        # checks that checkVarHGVS works for an equal nucleotide delins in BRCA2
        self.varBRCA2["Ref"] = "GC"
        self.varBRCA2["Alt"] = "AG"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.-26_-25delGCinsAG"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "TA"
        self.varBRCA2["Alt"] = "GC"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)

        
    def test_checkVarHGVSDelinsUnequal(self):
        '''Tests checkVarHGVS for delins variants with unequal number of nucleotides deleted and inserted'''
        
        # checks that checkVarHGVS works for an unequal nucleotide delins in BRCA1
        self.varBRCA1["Ref"] = "AGAAA"
        self.varBRCA1["Alt"] = "TT"
        self.varBRCA1["pyhgvs_cDNA"] = "NM_007294.3:c.-19-17_-19-13delTTTCTinsAA"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertTrue(checkVar)

        self.varBRCA1["Ref"] = "TCGAT"
        self.varBRCA1["Alt"] = "TA"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA1)
        self.assertNotEqual(True, checkVar)
        
        # checks that checkVarHGVS works for an unequal nucleotide delins in BRCA2
        self.varBRCA2["Ref"] = "TACCCCTATTG"
        self.varBRCA2["Alt"] = "ACAT"
        self.varBRCA2["pyhgvs_cDNA"] = "NM_000059.3:c.1232_1242delTACCCCTATTGinsACAT"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertTrue(checkVar)

        self.varBRCA2["Ref"] = "GCATGCTATGA"
        self.varBRCA2["Alt"] = "GTCA"
        checkVar = checkVarHGVS.checkVarHGVS(self.varBRCA2)
        self.assertNotEqual(True, checkVar)
