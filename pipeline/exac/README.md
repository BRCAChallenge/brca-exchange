### Build instructions for ExAC

1. Create a build directory, set environment variable $EXAC
2. pushd $EXAC; 
3. wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 
4. wget ftp://ftp.broadinstitute.org/pub/ExAC_release/current/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz.tbi
5. popd
6. tabix -h $EXAC/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 17:41191488-41322420 > $EXAC/exac.brca1.hg19.vcf
7. tabix -h $EXAC/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz 13:32889080-32973809 > $EXAC/exac.brca2.hg19.vcf
8. vcf-concat $EXAC/exac.brca1.hg19.vcf $EXAC/exac.brca2.hg19.vcf > $EXAC/exac.brca12.hg19.vcf
9. CrossMap.py vcf $BRCA_RESOURCES/hg19ToHg38.over.chain.gz $EXAC/exac.brca12.hg19.vcf $BRCA_RESOURCES/hg38.fa $EXAC/exac.brca12.hg38.vcf
10. vcf-sort $EXAC/exac.brca12.hg38.vcf > $EXAC/exac.brca12.sorted.hg38.vcf
