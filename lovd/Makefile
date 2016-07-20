releaseDir = $(BRCA_PIPELINE_DATA)
lovdDir = $(releaseDir)/data/lovd
vcfDir = $(releaseDir)/data/pipeline_input
leiden = ../leidenv1.0_package/

export PATH := .:$(leiden)/bin:$(PATH)
export PYTHONPATH := .:$(leiden):$(PYTHONPATH)


$(lovdDir)/sharedLOVD_brca1.vcf:	$(lovdDir)/BRCA1.txt
	./lovd2vcf -i $(lovDir)/BRCA1.txt -o $(lovdDir)/sharedLOVD_brca1.vcf -a sharedLOVDAnnotation -b 1

$(lovdDir)/sharedLOVD_brca2.vcf:	$(lovdDir)/BRCA2.txt
	./lovd2vcf -i BRCA2.txt -o sharedLOVD_brca2.vcf -a sharedLOVDAnnotation -b 1

$(lovdDir)/BRCA1.txt $(lovdDir)/BRCA2.txt:	Makefile
	python $(leiden)/bin/extract_data.py -u http://databases.lovd.nl/shared/ -l BRCA1 BRCA2 -o $(lovdDir)
