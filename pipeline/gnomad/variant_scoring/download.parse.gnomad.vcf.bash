#!/usr/bin/env bash

export BRCA1_CHROM=chr17
export BRCA2_CHROM=chr13

export BRCA1_START_HG19=41196312
export BRCA1_END_HG19=41277500
export BRCA2_START_HG19=32889617
export BRCA2_START_HG19=32973809

export BRCA1_START_HG38=43044295
export BRCA1_END_HG38=43125483
export BRCA2_START_HG38=32315474
export BRCA2_END_HG38=32400266


wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr17.vcf.bgz.tbi
bcftools view -r ${BRCA1_CHROM}:${BRCA1_START_HG38}-${BRCA1_END_HG38} gnomad.exomes.v4.1.sites.chr17.vcf.bgz \
	 -o brca1.gnomAD.exomes.4.1.hg38.vcf -Ov

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr13.vcf.bgz.tbi
bcftools view -r ${BRCA2_CHROM}:${BRCA2_START_HG38}-${BRCA2_END_HG38} gnomad.exomes.v4.1.sites.13.chr13.vcf.bgz \
	 -o brca2.gnomAD.exomes.4.1.hg38.vcf -Ov

vcf-concat brca1.gnomAD.exomes.4.1.hg38.vcf brca2.gnomAD.exomes.4.1.hg38.vcf \
	   > brca.gnomAD.exomes.4.1.hg38.vcf
rm *bgz*

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr17.vcf.bgz.tbi
bcftools view -r ${BRCA1_CHROM}:${BRCA1_START_HG38}-${BRCA1_END_HG38} gnomad.joint.v4.1.sites.chr17.vcf.bgz \
	 -o brca1.gnomAD.joint.4.1.hg38.vcf -Ov

wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr13.vcf.bgz.tbi
bcftools view -r ${BRCA2_CHROM}:${BRCA2_START_HG38}-${BRCA2_END_HG38} gnomad.joint.v4.1.sites.chr13.vcf.bgz \
	 -o brca2.gnomAD.joint.4.1.hg38.vcf -Ov

vcf-concat brca1.gnomAD.joint.4.1.hg38.vcf brca2.gnomAD.joint.4.1.hg38.vcf \
	   > brca.gnomAD.joint.4.1.hg38.vcf
rm *bgz*

extract_from_vcf.py -i brca.gnomAD.joint.4.1.hg38.vcf > brca.gnomAD.joint.4.1.hg38.for_Johanna.tsv


