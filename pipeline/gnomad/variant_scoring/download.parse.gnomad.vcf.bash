#!/usr/bin/env bash

export BRCA1_CHROM=17
export BRCA2_CHROM=13

export BRCA1_START_HG19=41196312
export BRCA1_END_HG19=41277500
export BRCA2_START_HG19=32889617
export BRCA2_START_HG19=32973809

export BRCA1_START_HG38=43044295
export BRCA1_END_HG38=43125483
export BRCA2_START_HG38=32315474
export BRCA2_END_HG38=32400266


wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.17.vcf.bgz.tbi
bcftools view -r ${BRCA1_CHROM}:${BRCA1_START_HG19}-${BRCA1_END_HG19} gnomad.exomes.r2.1.1.sites.17.vcf.bgz \
	 -o brca1.gnomAD.exomes.2.1.1.hg19.vcf -Ov

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.17.vcf.bgz.tbi
bcftools view -r ${BRCA1_CHROM}:${BRCA1_START_HG19}-${BRCA1_END_HG19} gnomad.genomes.r2.1.1.sites.17.vcf.bgz \
	 -o brca1.gnomAD.genomes.2.1.1.hg19.vcf -Ov

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.13.vcf.bgz.tbi
bcftools view -r ${BRCA2_CHROM}:${BRCA2_START_HG19}-${BRCA2_END_HG19} gnomad.exomes.r2.1.1.sites.13.vcf.bgz \
	 -o brca2.gnomAD.exomes.2.1.1.hg19.vcf -Ov

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.13.vcf.bgz.tbi
bcftools view -r ${BRCA2_CHROM}:${BRCA2_START_HG19}-${BRCA2_END_HG19} gnomad.genomes.r2.1.1.sites.13.vcf.bgz \
	 -o brca2.gnomAD.genomes.2.1.1.hg19.vcf -Ov

vcf-concat brca1.gnomAD.exomes.2.1.1.hg19.vcf brca2.gnomAD.exomes.2.1.1.hg19.vcf brca1.gnomAD.genomes.2.1.1.hg19.vcf \
	   brca2.gnomAD.genomes.2.1.1.hg19.vcf > brca.gnomAD.2.1.1.hg19.vcf
extract_from_vcf.py -i brca.gnomAD.2.1.1.hg19.vcf > brca.gnomAD.2.1.1.hg19.flags.tsv
rm *bgz*


wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr17.vcf.bgz.tbi
bcftools view -r chr${BRCA1_CHROM}:${BRCA1_START_HG38}-${BRCA1_END_HG38} gnomad.genomes.v3.1.1.sites.chr17.vcf.bgz \
	 -Ov -o brca1.gnomAD.genomes.3.1.1.hg38.vcf 
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr13.vcf.bgz.tbi
bcftools view -r chr${BRCA2_CHROM}:${BRCA2_START_HG38}-${BRCA2_END_HG38} gnomad.genomes.v3.1.1.sites.chr13.vcf.bgz \
	 -Ov -o brca2.gnomAD.genomes.3.1.1.hg38.vcf 

vcf-concat brca1.gnomAD.genomes.3.1.1.hg38.vcf brca2.gnomAD.genomes.3.1.1.hg38.vcf  > brca.gnomAD.3.1.1.hg38.vcf
extract_from_vcf.py -i brca.gnomAD.3.1.1.hg38.vcf > brca.gnomAD.3.1.1.hg38.flags.tsv
rm *bgz*

