'''
NOTE:
This file has several issues, and has therefore been renamed to remove it from testing.
To replace, change the name from testloadingscript.py to test_loadingscript.py.
'''

import nose2
import unittest
import logging
import pyhgvs
import pyhgvs.utils as pyhgvs_utils
from pygr.seqdb import SequenceFileDB
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import csv

# function for getting variants from built.tsv
def load_variants(variant_file):
    variants_tsv = open(variant_file)
    reader = csv.reader(variants_tsv, dialect="excel-tab")
    header = next(reader)
    variant_row = []
    
    for row in reader:
        row_dict = dict(list(zip(header, row)))
        variant_row.append(row_dict)

    return variant_row
    variants_tsv.close()


class VariantTSVTest(unittest.TestCase):


	def setUp(self):
		self.data = load_variants("built.tsv")

	#def tearDown(self):
	
	def test_data_not_empty(self):
		for i in self.data:
			for v in i.values():
				self.assertIsNot(v, "") 


	def test_data_isDict(self):
		for i in self.data:
			self.assertTrue(i)


	def test_chr_is_digit(self):
		for i in self.data:
			self.assertTrue(i['Chr'].isdigit()) #or .isalpha())


	def test_chr_is_13_or_17(self):
		chrom = ['13', '17']
		for i in self.data:
			self.assertIn(i['Chr'], chrom) 


	def test_gene_is_BRCA1_or_BRCA2(self):
		gene = ['BRCA1', 'BRCA2']
		for i in self.data:
			self.assertIn(i['Gene_Symbol'], gene)


	def test_pos_is_alphanumeric(self):
		for i in self.data:
			self.assertTrue(i['Pos'].isalnum()) 


	def test_ref_valid_nucleotide(self):
		nucleotide = ['A', 'T', 'C', 'G']
		for i in self.data:
			for n in i['Ref']: 
				self.assertIn(n, nucleotide) 

	# Alt maybe 'R' or 'Y' as an 'N'
	def test_alt_valid_nucleotide(self):
		nucleotide = ['A', 'T', 'C', 'G', 'R', 'Y']
		for i in self.data:
			for n in i['Alt']: 
				self.assertIn(n, nucleotide) 


	def test_pos_hgvs38start_are_consistent(self):
		for i in self.data:
			self.assertEqual(i['Pos'], i['Hg38_Start']) 


	def test_pos_genomic_coord_are_consistent(self):
		for i in self.data:
			gc = (i['Genomic_Coordinate_hg38'].split(":")[1].split(".")[1])
			if len(gc) == 10:
				genomic_coord = int(gc.split(".")[1])
			else:
				genomic_coord = int(gc)
			self.assertEqual(i['Pos'], str(genomic_coord))

	
	def test_pos_pyhgvscoord_are_consistent(self):
		for i in self.data:
			pyhgvs = (i['pyhgvs_Genomic_Coordinate_38'].split(":")[1].split(".")[1])
			self.assertEqual(i['Pos'], pyhgvs)
				

	def test_pyhgvs_cdna_coordinate_correct(self):
		for i in self.data:
			pyhgvs_coord = i['pyhgvs_Genomic_Coordinate_38']
			pyhgvs_cDNA = i['pyhgvs_cDNA']
			genome = SequenceFileDB('../reference_genome/hg38/hg38.fa')
			
			def get_transcript(name):
			    REFGENE = "../refgene38_brca.txt"
			    with open(REFGENE) as infile:
			        TRANSCRIPTS = pyhgvs_utils.read_transcripts(infile)
			    return TRANSCRIPTS.get(name)

			chrom, offset, ref, alt = pyhgvs.parse_hgvs_name(pyhgvs_cDNA, genome, 
				get_transcript=get_transcript)
			test_coord = chrom + ":" + "g." + str(offset) + ":" + ref + ">" + alt

			self.assertEqual(pyhgvs_coord, test_coord)

	# Errors in this test if 'R' or 'Y' in Alt allele
	def test_pyhgvs_cdna_protein_correct(self):
		HDP = hgvs.dataproviders.uta.connect()
		EVM = hgvs.variantmapper.EasyVariantMapper(HDP,
    			primary_assembly='GRCh37', alt_aln_method='splign')
		HP = hgvs.parser.Parser()
		for i in self.data:
			pyhgvs_cDNA = i['pyhgvs_cDNA']
			pyhgvs_Protein = i['pyhgvs_Protein'].split(":")[1]
			hgvs_c = HP.parse_hgvs_variant(pyhgvs_cDNA)
			hgvs_p = str(EVM.c_to_p(hgvs_c)).split(":")[1]
			self.assertEqual(pyhgvs_Protein, hgvs_p)


if __name__ == '__main__':
	unittest.main()

