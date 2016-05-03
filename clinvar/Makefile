releaseDir = /hive/groups/cgl/brca/release1.0
clinVarDir = $(releaseDir)/data/ClinVar
vcfDir = $(releaseDir)/pipeline_input

#export PATH := .:../data_merging:$PATH

$(vcfDir)/ClinVarBrca.vcf:	$(clinVarDir)/ClinVarBrca.vcf
	cp $< $@

$(clinVarDir)/ClinVarBrca.vcf:	$(clinVarDir)/ClinVarBrca.txt ../data_merging/convert_tsv_to_vcf.py
	../data_merging/convert_tsv_to_vcf.py -i $< -o $@ -s "ClinVar"

$(clinVarDir)/ClinVarBrca.txt:	$(clinVarDir)/ClinVarBrca.xml clinVarParse.py
	clinVarParse.py	$< --assembly GRCh38 > $@

$(clinVarDir)/ClinVarBrca.xml:	$(clinVarDir)/ClinVarFullRelease.xml clinVarBrca.py
	clinVarBrca.py $< |awk '{ if (NF > 0) { print }}' | > $@
