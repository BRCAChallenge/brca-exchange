#!/bin/bash

BASE=$1
IMAGE=$2
PREVIOUS_RELEASE=$3
RELEASE_NOTES=$4
DATA_DATE=$5

# optional: command to execute in container
RUN_COMMAND=$6

SYNAPSE_CACHE=$(mktemp -d /tmp/synapse_cache_pipeline_XXXXX)

trap "rm -r ${SYNAPSE_CACHE}" EXIT

docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) \
       -e "DATA_DATE=${DATA_DATE}" \
       -v ${BASE}/resources:/files/resources \
       -v ${BASE}/brca_out:/files/data \
       -v /home/pipeline/monthly_releases/scripts/luigi_pipeline_credentials.cfg:/opt/luigi_pipeline_credentials.cfg \
       -v ${PREVIOUS_RELEASE}:/files/previous_release.tar.gz \
       -v ${RELEASE_NOTES}:/files/release_notes.txt \
       -v ${SYNAPSE_CACHE}:/.synapseCache \
       ${IMAGE} ${RUN_COMMAND}
