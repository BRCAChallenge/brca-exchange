#!/usr/bin/env bash

if [ "$#" -ne "1" ]; then
    echo "Downloads required resources files such as reference sequences"
    echo "Usage: download_resource_files.sh [directory]"
    exit 2
fi

echo "Downloading into $1 ...."

cd $1

URL_PREFIX="https://brcaexchange.org/backend/downloads/resources"

HASH=shasum

declare -A files

files[brca1_hg19.txt]=0103381d9f931367e85cc7d4c0504973a6ef08e2
files[brca1_hg38.txt]=1adac34fa794b8057926b836384c52fe53934b6a
files[brca2_hg19.txt]=a0f76d6a61d5275be4695e9628b0bf8d49dc1a37
files[brca2_hg38.txt]=e3d40430167100bb3252a5f0e3655f56e244e413
files[hg18.fa.gz]=177dcd8262eed25f565cd2cd10f87fe63a7e54cb
files[hg19.fa.gz]=c1a06a3122bb9b3747ceeefb902aecd3ea9b9008
files[hg19ToHg38.over.chain.gz]=4769aaa7a5a7837c3947604207da346b1dd19a8f
files[hg38ToHg19.over.chain.gz]=48d40f5aa131f22ff41fad441a085d6c464ef56d
files[hg38.fa.gz]=d329507d63d32405664053785b2bfdc02b3e6b28
files[refseq_annotation.hg18.gp]=21911f9cd02f00b82eeb6068b7ab737c4d89f896
files[refseq_annotation.hg19.gp]=e12092220926fe112da26ffdc8188b76aef83231
files[refseq_annotation.hg38.gp]=9fb29159e516d2f28c2ec0e41a338c85ac699528

for file in "${!files[@]}"; do
    hash=$(${HASH} $file | cut -f1 -d' ')
    if [ "$hash" != "${files[$file]}" ]; then
        echo "Downloading ${file}..."
        wget ${URL_PREFIX}/${file} -O ${file}
    fi
    hash=$(${HASH} $file | cut -f1 -d' ')
    if [ "$hash" != "${files[$file]}" ]; then
        echo "Hash didn't match for ${file}"
        exit -1
    fi
done

for f in hg18.fa.gz hg19.fa.gz hg38.fa.gz; do
    echo "Unzipping ${f} ..."
    gunzip < "$f" > "$(basename $f '.gz')"
done
