releaseDir = $(BRCA_PIPELINE_DATA)
resources = $(releaseDir)/pipeline-resources
inputVcfDir = $(releaseDir)/pipeline-input
inputVcfs = $(inputVcfDir)/1000G_brca.sorted.hg38.vcf \
	$(inputVcfDir)/ClinVarBrca.vcf \
	$(inputVcfDir)/sharedLOVD_brca12.sorted.hg38.vcf \
	$(inputVcfDir)/exLOVD_brca12.sorted.hg38.vcf \
	$(inputVcfDir)/bic_brca12.sorted.hg38.vcf \
	$(inputVcfDir)/exac.brca12.sorted.hg38.vcf \
	$(inputVcfDir)/esp.brca12.sorted.hg38.vcf
outputFiles = $(releaseDir)/pipeline-output


all:	$(outputFiles)/aggregated.tsv

$(outputFiles)/aggregated.tsv:	$(outputFiles)/merged.withVEP.csv aggregate_across_columns.py brca_pseudonym_generator.py
	aggregate_across_columns.py -i $< -o $(outputFiles)/tmp.tsv
	brca_pseudonym_generator.py -i $(outputFiles)/tmp.tsv \
	  -j $(resources)/hg18.fa \
	  -k $(resources)/hg19.fa -l $(resources)/hg38.fa \
	  -r $(resources)/refseq_annotation.hg18.gp \
	  -s $(resources)/refseq_annotation.hg19.gp \
	  -t $(resources)/refseq_annotation.hg38.gp \
	  -o $@

$(outputFiles)/merged.withVEP.csv:	$(outputFiles)/merged.csv add_VEP_result.py
	add_VEP_result.py -i $< -v $(inputVcfDir)/vep_output_3_3_2016.vcf -o $@

$(outputFiles)/merged.csv:	$(inputVcfs) variant-merging.py
	variant-merging.py -i $(inputVcfDir)/ -o $(releaseDir)/pipeline_output/ -p


