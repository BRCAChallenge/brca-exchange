# download from http://research.nhgri.nih.gov/bic/
# username/password is in /hive/groups/cgl/brca/phase1/data/bic/account.txt at UCSC
# brca1 data
# https://research.nhgri.nih.gov/projects/bic/Member/brca1_mutation_database.shtml
# brca2 data 
# https://research.nhgri.nih.gov/projects/bic/Member/brca2_mutation_database.shtml
# datafiles are called brcaX_data.txt

# convert data to vcf, CURRENTLY ONLY EXPORTS THE SNPS !!
python convBic.py  > bicSnp.vcf
# index
bgzip bicSnp.vcf
tabix -p vcf bicSnp.vcf.gz 

# lift to hg38
CrossMap.py vcf /gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz bicSnp.vcf /hive/data/genomes/hg19/hg19.fa bicSnp.hg38.vcf
#@ 2015-05-13 09:28:15: Read chain_file:
#/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
#@ 2015-05-13 09:28:16: Load reference genome: /hive/data/genomes/hg19/hg19.fa
#@ 2015-05-13 09:31:06: Total entries: 19357
#@ 2015-05-13 09:31:06: Failed to map: 5253
#
# !! failed to map 5k variants!

# index
bgzip bicSnp.hg38.vcf
tabix -p vcf bicSnp.hg38.vcf.gz 
