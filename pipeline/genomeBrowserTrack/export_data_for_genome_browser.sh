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

HG19_BED=${TMP_OUT}/brcaExchange.hg19.bed
HG38_BED=${TMP_OUT}/brcaExchange.hg38.bed
AS=${TMP_OUT}/brcaExchange.as

python ${SRC_DIR}/brcaToBed.py -i ${TMP_DIR}/output/release/built_with_change_types.tsv -o19 ${HG19_BED} -o38 ${HG38_BED} -a ${AS}

sort -k1,1 -k2,2 ${HG19_BED} -o ${HG19_BED}
sort -k1,1 -k2,2 ${HG38_BED} -o ${HG38_BED}

set +e

# bedToBigBed errors on variants > 255 characters, this removes those problem variants until bedToBigBed succeeds

# captures error message if present
ERROR=`bedToBigBed -type=bed9+ -as=${AS} -tab ${HG19_BED} ${SRC_DIR}/hg19.chrom.sizes ${TMP_HG19}/brcaExchange.bb 2>&1 `

# parses line number of problem variant from error
ERRORLINE=`echo $ERROR | sed -n -e 's/^.*line //p' | sed 's/\s.*$//'`

until [ -z "$ERRORLINE" ]
do
	# removes error variant from both BED files if present
	sed -i "${ERRORLINE}d" ${HG19_BED}
	sed -i "${ERRORLINE}d" ${HG38_BED}
	echo "removed variant due to error:"
	echo $ERROR

	# captures error message if present
        ERROR=`bedToBigBed -type=bed9+ -as=${AS} -tab ${HG19_BED} ${SRC_DIR}/hg19.chrom.sizes ${TMP_HG19}/brcaExchange.bb 2>&1 `
        echo $ERROR

        # parses line number of problem variant from error
        ERRORLINE=`echo $ERROR | sed -n -e 's/^.*line //p' | sed 's/\s.*$//'`
        echo $ERRORLINE
done

set -e

bedToBigBed -type=bed9+ -as=${AS} -tab ${HG19_BED} ${SRC_DIR}/hg19.chrom.sizes ${TMP_HG19}/brcaExchange.bb
bedToBigBed -type=bed9+ -as=${AS} -tab ${HG38_BED} ${SRC_DIR}/hg38.chrom.sizes ${TMP_HG38}/brcaExchange.bb

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
