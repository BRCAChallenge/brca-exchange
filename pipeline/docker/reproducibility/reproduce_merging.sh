#!/bin/bash

###############################################################################
# Script to partially reproduce the result from the BRCA Exchange pipeline.
###############################################################################

if [[ $# -ne 3 ]]; then
    echo "Expecting 3 arguments, got $#"
    echo "Usage $0: base_dir release_date previous_release_date"
    echo ""
    echo "where:"
    echo "    base_dir: path to base directory for the analysis"
    echo "    release_date: date of the release to reproduce, e.g. '2018-02-17'"
    echo "    previous_release_date: date of the release prior to release to be reproduced, e.g. '2018-01-16'"
    echo ""
    echo "Example:"
    echo ""
    echo "./reproduce_merging.sh /tmp/reproductions/ 2018-02-17 2018-01-16"
    echo ""
    echo "Note: you can lookup release dates here: https://brcaexchange.org/releases"
    exit 2
fi

BASE=$1
RELEASE_DATE=$2
PREVIOUS_RELEASE_DATE=$3

IMAGE="brcachallenge/brca-exchange-pipeline:data_release_${RELEASE_DATE}"

###############################################################################

# convert to absolute path (required by docker daemon)
if [[ ! "${BASE}" =~ ^/ ]]
then
    BASE="$(pwd)/${BASE}"
fi

SED_DATE_CONVERSION_CMD="s/20([0-9]{2})-([0-9]{2})-([0-9]{2})/\2-\3-\1/g"

function download_release_archive()
{
    DATE_US=$1
    TARGET=$2

    URL="https://brcaexchange.org/backend/downloads/releases/release-${DATE_US}/release-${DATE_US}.tar.gz"
    wget "${URL}" -O ${TARGET}

    if [ "$?" -ne "0" ]; then
        echo "Failed to download archive $URL. Please check https://brcaexchange.org/releases for the available release dates"
        exit 1
    fi
}

RELEASE_DATE_US=$(echo ${RELEASE_DATE} | sed -E ${SED_DATE_CONVERSION_CMD})
RELEASE_ARCHIVE=${BASE}/release-${RELEASE_DATE_US}.tar.gz
download_release_archive ${RELEASE_DATE_US} ${RELEASE_ARCHIVE}

PREVIOUS_RELEASE_DATE_US=$(echo ${PREVIOUS_RELEASE_DATE} | sed -E ${SED_DATE_CONVERSION_CMD})
PREVIOUS_RELEASE_ARCHIVE=${BASE}/release-${PREVIOUS_RELEASE_DATE_US}.tar.gz
download_release_archive ${PREVIOUS_RELEASE_DATE_US} ${PREVIOUS_RELEASE_ARCHIVE}

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

OUTPUT_DIR=${BASE}/data_out
mkdir -p ${OUTPUT_DIR}

# extracting required files from release tar
# excluding the output/release directory, because this is what we want to reproduce
[[ "$OSTYPE" == "linux-gnu" ]] && WILDCARDS_OPT='--wildcards' # hacking around differing support across platforms of the --wildcards tar option
tar zxf ${RELEASE_ARCHIVE} ${WILDCARDS_OPT} -C ${OUTPUT_DIR} --exclude 'output/release' 'output/*.vcf*' 'output/*.tsv'

echo "$(date) Running BRCA Exchange pipeline (merging part)"
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) \
       -e "DATA_DATE=${RELEASE_DATE}" \
       -v ${RESOURCE_DIR}:/files/resources \
       -v ${OUTPUT_DIR}:/files/data \
       -v ${PREVIOUS_RELEASE_ARCHIVE}:/files/previous_release.tar.gz \
       -v ${RELEASE_NOTES}:/files/release_notes.txt \
       ${IMAGE}

echo "If everything went well, find the reproduced archive in ${OUTPUT_DIR}/release-${RELEASE_DATE_US}.tar.gz"
