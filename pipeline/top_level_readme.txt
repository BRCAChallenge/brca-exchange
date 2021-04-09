(Readme last updated April 06, 2021)

BRCA Exchange download data

This archive contains the variant data released by BRCA Exchange
(brcaexchange.org) for one particular release, along with the pre-processed
data from which it was aggregated.  This file describes the contents of that
archive.  If you have any questions or concerns, please feel free to contact
brca-exchange-contact@genomicsandhealth.org.

A good starting point to explore the data is the `output/variants_output.tsv` file. Have a look at the Jupyter notebook in
https://github.com/BRCAChallenge/brca-exchange/tree/master/pipeline/notebooks/variants_output_file_tour.ipynb
to see how the data can be used. The columns in the file are described in the `output/variants_output_field_metadata.tsv` file.

Description of the files in the tar ball:

output/
	Top-level directory, containing pre-processed VCFs of variant data
	from the contributing repositories, filtered for relevant variants,
	and translated into GRCh38 coordinates as needed.
output/1000G.sorted.hg38.vcf
	Variants selected from 1000 Genomes, with genotype data.
output/ClinVar.vcf
	Variants from ClinVar at NCBI.
output/bic_brca12.sorted.hg38.vcf
	Variants from the Breast Cancer Information Core (BIC) at NCBI.
output/enigma_from_clinvar.tsv
	Tab-delimited file of the ENIGMA clinical significance data.
output/esp.sorted.hg38.vcf
	Variants from the NCBI Exome Sequencing Project
output/exLOVD_brca12.sorted.hg38.vcf
	Variants from the BRCA ExUV project at the Huntsman Institute
output/exac.brca12.sorted.hg38.vcf
	Variants from the ExAC project at the Broad Institute
output/findlay_BRCA1_ring_function_scores.clean.sorted.hg38.vcf
	Functional assay scores from Findlay et al, 2018.
output/md5sums.txt
	md5sums for all data files in this archive
output/release/
	Subdirectory containing the computed data released on the BRCA Exchange web portal
output/release/artifacts/
	Subdirectory with miscellaneous intermediate results from the data aggregation pipeline
output/release/built_final.tsv
	Final and complete output of pipeline. Currently just a symlink to output/release/built_with_change_types.tsv
output/release/built_with_change_types.tsv
	Tab-delimited file of the variants shared on BRCA Exchange, their attributes, and their update status since the
	last release.
output/release/diff/
	Subdirectory containing files that detail how the set of variants has changed since the previous release
output/release/diff/added.tsv
	Tab-delimited file containing new variants
output/release/diff/added_data.tsv
	Tab-delimited file containing updated variants
output/release/diff/diff.json
	JSON file detailing specific changes to the updated variants
output/release/diff/diff.txt
	Text file detailing specific changes to the updated variants, with old and new data values
output/release/diff/removed.tsv
	Tab-delimited file containing any variants removed in the last release
output/release/field_metadata.tsv
	File documenting the columns of built_final.tsv
output/release/metadata/
	Subdirectory containing metadata on this release
output/release/metadata/version.json
	JSON file containing metadata on this release.
output/sharedLOVD.sorted.hg38.vcf
	Variants from the shared LOVD repository at Leiden University
output/variants_output.tsv
    Simplified version of built_final.tsv for downstream consumption
output/variants_output_field_metadata.tsv
    File documenting the columns of variants_output.tsv

For any concerns with this data, please contact
brca-exchange-contact@genomicsandhealth.org.

