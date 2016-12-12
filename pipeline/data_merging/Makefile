releaseDir = $(BRCA_PIPELINE_DATA)
resources = $(releaseDir)/pipeline-resources/
inputVcfDir = $(releaseDir)/pipeline-input
inputVcfs = $(inputVcfDir)/1000G_brca.sorted.hg38.vcf \
	$(inputVcfDir)/ClinVarBrca.vcf \
	$(inputVcfDir)/sharedLOVD_brca12.sorted.hg38.vcf \
	$(inputVcfDir)/exLOVD_brca12.sorted.hg38.vcf \
	$(inputVcfDir)/bic_brca12.sorted.hg38.vcf \
	$(inputVcfDir)/exac.brca12.sorted.hg38.vcf \
	$(inputVcfDir)/esp.brca12.sorted.hg38.vcf
outputFiles = $(releaseDir)/pipeline-output/


all:	$(outputFiles)/built.tsv

$(outputFiles)/built.tsv:	$(outputFiles)/aggregated.tsv brca_pseudonym_generator.py
	brca_pseudonym_generator.py -i $< -p \
	  -j $(resources)/hg18.fa \
	  -k $(resources)/hg19.fa -l $(resources)/hg38.fa \
	  -r $(resources)/refseq_annotation.hg18.gp \
	  -s $(resources)/refseq_annotation.hg19.gp \
	  -t $(resources)/refseq_annotation.hg38.gp \
	  -o $@

$(outputFiles)/aggregated.tsv:	$(outputFiles)/annotated.tsv aggregate_across_columns.py 
	aggregate_across_columns.py -i $< -o $@

$(outputFiles)/annotated.tsv:	$(outputFiles)/merged.tsv add_annotation.py
	add_annotation.py -i $< -o $@

$(outputFiles)/merged.tsv:	$(inputVcfs) variant_merging.py
	variant_merging.py -i $(inputVcfDir)/ -o $(outputFiles) -p -r $(resources)


