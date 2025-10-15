#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "need one argument: data destination directory"
    exit 1
fi

DEST_DIR=$1

if [ ! -d ${DEST_DIR} ]; then
    echo "create ${DEST_DIR}"
    mkdir -p ${DEST_DIR}
fi

echo "Start downloading into ${DEST_DIR}"

cd ${DEST_DIR}

echo "Download vcf for chr 13"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.13.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr13.vcf.bgz


echo "Download vcf for chr 17"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.17.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.chr17.vcf.bgz

echo "Download coverage summaries"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.summary.tsv.bgz



