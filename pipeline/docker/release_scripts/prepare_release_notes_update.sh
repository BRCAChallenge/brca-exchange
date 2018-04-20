#!/bin/bash

if [ "$#" -ne "1" ]; then
    echo "Expecting release tag, i.e. data_release_yyyyMMdd, as argument!"
    exit 1
fi

RELEASE_TAG=$1

DATA_DIR="/home/pipeline/monthly_releases/${RELEASE_TAG}/brca_out"

# deleting some files, s.t. when running the luigi pipeline again, the updated release notes will be included.
rm ${DATA_DIR}/release-*.tar.gz
rm ${DATA_DIR}/output/md5sums.txt
rm ${DATA_DIR}/output/README.txt
rm ${DATA_DIR}/output/release/metadata/version.json

