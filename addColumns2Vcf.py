"""
this script add three more INFO fields to a VCF file
these three INFO fields are:
HGVS_G: HGVS genomic coordinate
HGVS_C: HGVS cDNA coordinate(s)
HGVS_P: HGVS protein coordinate(s)

Note: this script can only handle SNPs, cann't do indels
"""

import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.variantmapper
import random



def main():
	old_vcf = open("data/brca.clinvar.10line.addHGVS.vcf", "r")
	new_vcf = open("data/brca.clinvar.mockdata.vcf", "w")
	new_info = ["FREQ", "BIC_P", "BIC_N", "MUTTYPE", "IARC", "DBSource"]


	for line in old_vcf:
		if "#" in line:
			new_vcf.write(line)
		else:
			line = add_info(line)
			new_vcf.write(line)

def add_info(line):
	current_info = line.strip().split("\t")[-1]
	current_info += ";FREQ=" + str(random.random())
	current_info += ";BIC_P=?;BIC_N=?"
	current_info += ";MUTTYPE=" + random.choice(["missense", "nonsense", "silent"])
	current_info += ";IARC=" + random.choice(["group 1", "group 2a", "group 2b", "group 3", "group 4"])
	current_info += ";DBSource=" + random.choice(["Clinvar", "LOVD", "BIC", "ExAC", "1000Genomes"])
	current_info += "\n"
	return current_info


def add_HGVS(line):
	items = line.strip().split("\t")
	chrom = items[0]
	pos = items[1]
	ref = items[3]
	alt = items[4]

	hgvs_g = "NC_0000" + str(chrom) + ".10:g." + str(pos) + ref + ">" + alt
	print hgvs_g
	hp = hgvs.parser.Parser()
	hgvs_g = hp.parse_hgvs_variant(hgvs_g)
	hdp = hgvs.dataproviders.uta.connect()
	evm = hgvs.variantmapper.EasyVariantMapper(hdp,
		primary_assembly='GRCh37', alt_aln_method='splign')
	transcripts = evm.relevant_transcripts(hgvs_g)
	for transcript in transcripts:
		hgvs_c = evm.g_to_c(hgvs_g, transcript)
		hgvs_p = evm.c_to_p(hgvs_c)
	
	items[-1] += ";HGVS_G=" + str(hgvs_g) + ";HGVS_C=" + str(hgvs_c) + ";HGVS_P=s" + str(hgvs_p)
	new_line = "\t".join(items) + "\n"
	return new_line

if __name__ == "__main__":
	main()








