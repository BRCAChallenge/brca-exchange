(Readme last updated December 15, 2017)

BRCA Exchange download data

This archive contains the variant data released by BRCA Exchange 
(brcaexchange.org) for one particular release, along with the pre-processed
data from which it was aggregated.  This file describes the contents of that
archive.

output/  Top-level directory, containing pre-processed VCFs of variant data
         from the contributing repositories, filtered for the BRCA variants,
	 and translated into GRCh38 coordinates as needed.
output/1000G_brca.sorted.hg38.vcf   		Variants selected from 1000  
				    		Genomes, with genotype data.
output/1000G_brca.sorted.hg38.vcffor_pipeline	Variants from 1000 Genomes,
						without genotype data
output/ClinVarBrca.vcf			BRCA variants from ClinVar at NCBI.
output/ENIGMA_combined_with_bx_ids.tsv	Tab-delimited file of the ENIGMA
					clinical significance data.   
output/bic_brca12.sorted.hg38.vcf	Variants from the Breast Cancer 
					Information Core (BIC) at NCBI.
output/esp.brca12.sorted.hg38.vcf	Variants from the NCBI Exome Sequencing
					Project
output/exLOVD_brca12.sorted.hg38.vcf	Variants from the BRCA ExUV project
					at the Huntsman Institute
output/exac.brca12.sorted.hg38.vcf	Variants from the ExAC project at the
					Broad Institute
output/md5sums.txt			md5sums for all data files in this 
					archive
output/release/				(see below)
output/sharedLOVD_brca12.sorted.hg38.vcf     BRCA variants from the shared
					     LOVD repository at Leiden 
					     University

output/release/	Subdirectory containing the computed data released on the
		BRCA Exchange web portal
output/release/artifacts/     Subdirectory with miscellaneous intermediate
			      results from the data aggregation pipeline
output/release/built.tsv      Tab-delimited file with the variants shared
			      on BRCA Exchange, and their attribute data
output/release/built_with_change_types.tsv	Same as above, with the 
						addition of a column that
						indicates if each variant is
						new or changed.

output/release/diff/	Subdirectory containing files that detail how the 
			set of variants has changed since the previous release
output/release/diff/added.tsv  		Tab-delimited file containing new 
					variants
output/release/diff/added_data.tsv	Tab-delimited file containing
					updated variants
output/release/diff/diff.json		JSON file detailing specific changes
					to the updated variants
output/release/diff/diff.txt		Text file detailing specific changes
					to the updated variants, with old and
					new data values
output/release/diff/removed.tsv		Tab-delimited file containing any
					variants removed in the last release


output/release/metadata/		Subdirectory containing metadata 
					on this	release
output/release/metadata/version.json	JSON file containing metadata on this
					release.

For any concerns with this data, please contact 
brca-exchange-contact@genomicsandhealth.org.






