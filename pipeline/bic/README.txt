# download from http://research.nhgri.nih.gov/bic/
# username/password is in /hive/groups/cgl/brca/phase1/data/bic/account.txt at UCSC
# brca1 data
# https://research.nhgri.nih.gov/projects/bic/Member/brca1_mutation_database.shtml
# brca2 data
# https://research.nhgri.nih.gov/projects/bic/Member/brca2_mutation_database.shtml
# datafiles are called brcaX_data.txt
# Put them in a working directory pointed to with the environment variable BIC.

# convert data to vcf, CURRENTLY ONLY EXPORTS THE SNPS !!
pushd $BRCA_RESOURCES; wget cline@hgwdev.soe.ucsc.edu:public_html/BRCA/resources/bicAnnotation; popd
./bic2vcf -i $BIC/brca1_data.txt -o $BIC/bic_brca1.hg19.vcf -b 1 -g $BRCA_RESOURCES/hg19.fa -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -a $BRCA_RESOURCES/bicAnnotation
./bic2vcf -i $BIC/brca2_data.txt -o $BIC/bic_brca2.hg19.vcf -b 2 -g $BRCA_RESOURCES/hg19.fa -r $BRCA_RESOURCES/refseq_annotation.hg19.gp -a $BRCA_RESOURCES/bicAnnotation
vcf-concat $BIC/bic_brca1.hg19.vcf $BIC/bic_brca2.hg19.vcf > $BIC/bic_brca12.hg19.vcf
CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $BIC/bic_brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $BIC/bic_brca12.hg38.vcf
vcf-sort $BIC/bic_brca12.hg38.vcf > $BIC/bic_brca12.sorted.hg38.vcf
cp $BIC/bic_brca12.sorted.hg38.vcf $PIPELINE_INPUT
