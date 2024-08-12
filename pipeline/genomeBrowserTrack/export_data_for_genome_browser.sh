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
TMP_HG19="${TMP_DIR}/hg19"
TMP_HG38="${TMP_DIR}/hg38" 
mkdir -p ${TMP_OUT}
mkdir -p ${TMP_HG19}
mkdir -p ${TMP_HG38}

SRC_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # path of this script plus some additional files
cp ${SRC_DIR}/hg19/trackDb.txt ${TMP_HG19}
cp ${SRC_DIR}/hg38/trackDb.txt ${TMP_HG38}

tar zxf ${RELEASE_ARCHIVE} -C ${TMP_DIR}

HG19_VAR_BED=${TMP_OUT}/variants.hg19.bed
HG38_VAR_BED=${TMP_OUT}/variants.hg38.bed
HG19_SV_BED=${TMP_OUT}/structural_variants.hg19.bed
HG38_SV_BED=${TMP_OUT}/structural_variants.hg38.bed
AS=${TMP_OUT}/brcaExchange.as

python ${SRC_DIR}/brcaToBed.py \
       -i ${TMP_DIR}/output/release/built_with_change_types.tsv \
       -l 50 \
       --output_hg19_var ${HG19_VAR_BED} \
       --output_hg38_var ${HG38_VAR_BED} \
       --output_hg19_sv ${HG19_SV_BED} \
       --output_hg38_sv ${HG38_SV_BED} \
       -a ${AS}


set +e

${SRC_DIR}/bigBedFromBed.sh ${HG19_VAR_BED} ${AS} ${SRC_DIR}/hg19.chrom_sizes ${TMP_HG19}/variants.bb
${SRC_DIR}/bigBedFromBed.sh ${HG38_VAR_BED} ${AS} ${SRC_DIR}/hg38.chrom_sizes ${TMP_HG38}/variants.bb
${SRC_DIR}/bigBedFromBed.sh ${HG19_SV_BED} ${AS} ${SRC_DIR}/hg19.chrom_sizes ${TMP_HG19}/structural_variants.bb
${SRC_DIR}/bigBedFromBed.sh ${HG38_SV_BED} ${AS} ${SRC_DIR}/hg38.chrom_sizes ${TMP_HG38}/structural_variants.bb

# preparing trackhubs directory, which then gets pushed to the remote machine
TRACKHUBS=${TMP_DIR}/trackhubs
[ ! -e ${TRACKHUBS} ] && mkdir ${TRACKHUBS}
cp -R ${SRC_DIR}/hub.txt ${SRC_DIR}/genomes.txt ${SRC_DIR}/description.html ${TRACKHUBS}
mkdir -p ${TRACKHUBS}/hg19
cp ${TMP_HG19}/* ${TRACKHUBS}/hg19
mkdir -p ${TRACKHUBS}/hg38
cp ${TMP_HG38}/* ${TRACKHUBS}/hg38



rsync -rtz --del --progress ${TRACKHUBS}/ brca@${HOST}:/var/www/html/production/trackhubs/

#rm -r ${TMP_DIR}
