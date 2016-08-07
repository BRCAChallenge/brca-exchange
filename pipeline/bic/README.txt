# download from http://research.nhgri.nih.gov/bic/
# username/password is in /hive/groups/cgl/brca/phase1/data/bic/account.txt at UCSC
# brca1 data
# https://research.nhgri.nih.gov/projects/bic/Member/brca1_mutation_database.shtml
# brca2 data 
# https://research.nhgri.nih.gov/projects/bic/Member/brca2_mutation_database.shtml
# datafiles are called brcaX_data.txt
# Put them in a working directory pointed to with the environment variable BIC.

# convert data to vcf, CURRENTLY ONLY EXPORTS THE SNPS !!
python convBic.py  --brca1 $BIC/brca1_data.txt --brca2 $BIC/brca2_data.txt > $BIC/bicSnp.hg19.vcf
# index
bgzip bicSnp.vcf
tabix -p vcf bicSnp.vcf.gz 

# lift to hg38
CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $BIC/bicSnp.hg19.vcf $BRCA_RESOURCES/hg38.fa $BIC/bicSnp.hg38.vcf
vcf-sort $BIC/bicSnp.hg38.vcf > $BIC/bicSnp.sorted.hg38.vcf
cp $BIC/bicSnp.sorted.hg38.vcf $PIPELINE_INPUT
