#!/usr/bin/env bash

# Generates bed files out of a BRCA Exchange release tar
#
# Usage: export_data_for_genome_browser.sh [release_tar]
#
# Note: requires bedToBigBed in PATH. bedToBigBed can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/ . Don't
# forget to `chmod +x` the file you download!

set -o errexit

HOST=${HOST:-brcaexchange.cloudapp.net}

RELEASE_ARCHIVE=$1
TMP_DIR="/tmp/export_genome_browser"
TMP_OUT="${TMP_DIR}/brcaToBed_out"

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # path of this script plus some additional files
cd ${SRC_DIR}

mkdir -p ${TMP_OUT}

tar zxf ${RELEASE_ARCHIVE} -C ${TMP_DIR}

HG19_BED=${TMP_OUT}/brcaExchange.hg19.bed
HG38_BED=${TMP_OUT}/brcaExchange.hg38.bed
AS=${TMP_OUT}/brcaExchange.as

python brcaToBed.py -i ${TMP_DIR}/output/release/built_with_change_types.tsv -o19 ${HG19_BED} -o38 ${HG38_BED} -a ${AS}

sort -k1,1 -k2,2 ${HG19_BED} -o ${HG19_BED}
sort -k1,1 -k2,2 ${HG38_BED} -o ${HG38_BED}

bedToBigBed -type=bed9+ -as=${AS} -tab ${HG19_BED} hg19.chrom.sizes hg19/brcaExchange.bb
bedToBigBed -type=bed9+ -as=${AS} -tab ${HG38_BED} hg38.chrom.sizes hg38/brcaExchange.bb

# preparing trackhubs directory, which then gets pushed to the remote machine
TRACKHUBS=${TMP_DIR}/trackhubs
[ ! -e ${TRACKHUBS} ] && mkdir ${TRACKHUBS}

cp -R hub.txt genomes.txt hg19 hg38 ${TRACKHUBS}

scp -rq ${TRACKHUBS}/. brca@${HOST}:/var/www/html/production/trackhubs/

rm -r ${TMP_DIR}
