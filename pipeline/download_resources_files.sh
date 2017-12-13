#!/bin/bash

if [ "$#" -ne "1" ]; then
    echo "Downloads required resources files such as reference sequences"
    echo "Usage: download_resource_files.sh [directory]"
    exit 2
fi

echo "Downloading into $1 ...."

cd $1


wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/refseq_annotation.hg38.gp
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/refseq_annotation.hg19.gp
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/refseq_annotation.hg18.gp
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/hg19.fa.gz && gunzip hg19.fa.gz
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/hg18.fa.gz && gunzip hg18.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/brca1_hg38.txt
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/brca2_hg38.txt
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/brca1_hg19.txt
wget http://hgwdev.soe.ucsc.edu/~cline/BRCA/resources/brca2_hg19.txt

