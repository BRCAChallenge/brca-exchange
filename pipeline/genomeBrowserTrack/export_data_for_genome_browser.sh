#!/usr/bin/env bash

# Generates bed files out of a BRCA Exchange release tar
#
# Usage: export_data_for_genome_browser.sh [release_tar]
#
# Note: requires bedToBigBed in PATH. bedToBigBed can be downloaded from http://hgdownload.cse.ucsc.edu/admin/exe/ . Don't
# forget to `chmod +x` the file you download!

set -o errexit

HOST=${HOST:-brcaexchange.org}

RELEASE_ARCHIVE=$1
TMP_DIR="/tmp/export_genome_browser"
TMP_OUT="${TMP_DIR}/brcaToBed_out"

mkdir -p ${TMP_OUT}

tar zxf ${RELEASE_ARCHIVE} -C ${TMP_DIR}

HG19_BED=${TMP_OUT}/brcaExchange.hg19.bed
HG38_BED=${TMP_OUT}/brcaExchange.hg38.bed
AS=${TMP_OUT}/brcaExchange.as

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # path of this script plus some additional files
python ${SRC_DIR}/brcaToBed.py -i ${TMP_DIR}/output/release/built_with_change_types.tsv -o19 ${HG19_BED} -o38 ${HG38_BED} -a ${AS}

sort -k1,1 -k2,2 ${HG19_BED} -o ${HG19_BED}
sort -k1,1 -k2,2 ${HG38_BED} -o ${HG38_BED}

bedToBigBed -type=bed9+ -as=${AS} -tab ${HG19_BED} ${SRC_DIR}/hg19.chrom.sizes ${SRC_DIR}/hg19/brcaExchange.bb
bedToBigBed -type=bed9+ -as=${AS} -tab ${HG38_BED} ${SRC_DIR}/hg38.chrom.sizes ${SRC_DIR}/hg38/brcaExchange.bb

# preparing trackhubs directory, which then gets pushed to the remote machine
TRACKHUBS=${TMP_DIR}/trackhubs
[ ! -e ${TRACKHUBS} ] && mkdir ${TRACKHUBS}

cp -R ${SRC_DIR}/hub.txt ${SRC_DIR}/genomes.txt ${SRC_DIR}/hg19 ${SRC_DIR}/hg38 ${TRACKHUBS}

rsync -rtz --del --progress ${TRACKHUBS}/ brca@${HOST}:/var/www/html/production/trackhubs/

rm -r ${TMP_DIR}
