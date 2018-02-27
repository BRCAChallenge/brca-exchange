#!/bin/bash

###############################################################################
# Script to partially reproduce the result from the BRCA Exchange pipeline.
#
# Please edit the PLEASE EDIT parts accordingly
###############################################################################

# Base directory
# PLEASE EDIT
BASE=~/please_edit

# date of the release to reproduce
# you can lookup release dates here: http://brcaexchange.org/releases
# PLEASE EDIT
RELEASE_DATE="2018-02-17" 

# previous release date of the release to reproduce
# PLEASE EDIT
PREVIOUS_RELEASE_DATE="2018-01-16" 

IMAGE="brcachallenge/brca-exchange-pipeline:data_release_${RELEASE_DATE}"

###############################################################################

SED_DATE_CONVERSION_CMD="s/20([0-9]{2})-([0-9]{2})-([0-9]{2})/\2-\3-\1/g"

RELEASE_DATE_US=$(echo ${RELEASE_DATE} | sed -E ${SED_DATE_CONVERSION_CMD})
RELEASE_ARCHIVE=${BASE}/release-${RELEASE_DATE_US}.tar.gz
wget "http://brcaexchange.org/backend/downloads/releases/release-${RELEASE_DATE_US}/release-${RELEASE_DATE_US}.tar.gz" -O ${RELEASE_ARCHIVE}

PREVIOUS_RELEASE_DATE_US=$(echo ${PREVIOUS_RELEASE_DATE} | sed -E ${SED_DATE_CONVERSION_CMD})
PREVIOUS_RELEASE_ARCHIVE=${BASE}/release-${PREVIOUS_RELEASE_DATE_US}.tar.gz
wget "http://brcaexchange.org/backend/downloads/releases/release-${PREVIOUS_RELEASE_DATE_US}/release-${PREVIOUS_RELEASE_DATE_US}.tar.gz" -O ${PREVIOUS_RELEASE_ARCHIVE}

# Download auxiliary resources
RESOURCE_DIR=${BASE}/resources
if [ ! -d ${RESOURCE_DIR} ]
then
    mkdir -p ${RESOURCE_DIR}

    echo "$(date) Downloading resources files..."
    docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) \
           -v ${RESOURCE_DIR}:/files/resources \
           ${IMAGE} \
           /opt/brca-exchange/pipeline/download_resources_files.sh /files/resources
fi

# Creating dummy release notes file
RELEASE_NOTES=${BASE}/some_release_notes.txt
touch ${RELEASE_NOTES}

# Creating dummy credentials file
CREDENTIALS_FILE=${BASE}/dummy_credentials_file.cfg
touch ${CREDENTIALS_FILE}

OUTPUT_DIR=${BASE}/brca_out
mkdir -p ${OUTPUT_DIR}

# extracting required files from release tar
# excluding the output/release directory, because this is what we want to reproduce
tar zxf ${RELEASE_ARCHIVE} -C ${OUTPUT_DIR} --exclude 'output/release' 'output/*.vcf*' 'output/*.tsv'

echo "$(date) Running BRCA Exchange pipeline (merging part)"
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) \
       -e "DATA_DATE=${RELEASE_DATE}" \
       -v ${RESOURCE_DIR}:/files/resources \
       -v ${OUTPUT_DIR}:/files/data \
       -v ${PREVIOUS_RELEASE_ARCHIVE}:/files/previous_release.tar.gz \
       -v ${RELEASE_NOTES}:/files/release_notes.txt \
       ${IMAGE}

echo "If everything went well, find the reproduced archive in ${OUTPUT_DIR}/release-${RELEASE_DATE_US}.tar.gz"
