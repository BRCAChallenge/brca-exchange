#!/usr/bin/env bash
#
# Generates an output BigBed file from a single input BED file,
# correcting for errors. bedToBigBed errors on variants > 255 characters.
# This script removes those problem variants until bedToBigBed succeeds
#
set -o errexit

INPUT_BED=$1
AS=$2
CHROM_SIZES=$3
OUTPUT_BIGBED=$4

sort -k1,1 -k2,2 ${INPUT_BED} -o ${INPUT_BED}

set +e

# captures the first error message if present
ERROR=`bedToBigBed -type=bed9+ -as=${AS} -tab ${INPUT_BED} ${CHROM_SIZES}  ${OUTPUT_BIGBED} 2>&1`

# parses line number of the first problem variant from the error message
ERRORLINE=`echo $ERROR | sed -n -e 's/^.*line //p' | sed 's/\s.*$//'`

until [ -z "$ERRORLINE" ]
do
    # removes error variant 
    sed -i "${ERRORLINE}d" ${INPUT_BED}
    echo $ERROR

    # captures and parses the line number of the next problem variant
    ERROR=`bedToBigBed -type=bed9+ -as=${AS} -tab ${INPUT_BED} ${CHROM_SIZES}  ${OUTPUT_BIGBED} 2>&1`
    ERRORLINE=`echo $ERROR | sed -n -e 's/^.*line //p' | sed 's/\s.*$//'`
    echo $ERRORLINE
done

set -e

bedToBigBed -type=bed9+ -as=${AS} -tab ${INPUT_BED} ${CHROM_SIZES} ${OUTPUT_BIGBED}
