#!/bin/bash

ARTIFACT_DIR="$1"
INPUT_FILE="$2"
OUTPUT_FILE="$3"
SEQ_REPO_DIR="${4:-/usr/local/share/seqrepo}"


# --- pre-step: build the seqrepo-rest-service docker
# See https://github.com/biocommons/seqrepo-rest-service

# ----------------------------------------------------
# --- 1. spin up seqrepo-rest-service, if needed
# ----------------------------------------------------

[ `docker ps -f name="seqrepo-rest-service" | wc -l` -gt 1 ] \
    || ../utilities/launch_seqrepo.sh ${SEQ_REPO_DIR}


# ...and wait for it to be available
./wait-for localhost:5000 && echo " => ...seqrepo-rest-service is ready"


# ----------------------------------------------------
# --- 2. execute append-vr-ids
# ----------------------------------------------------
PATH=../utilities:${PATH} python3 appendVRIds.py \
    -i ${ARTIFACT_DIR}/${INPUT_FILE} \
    -o ${ARTIFACT_DIR}/${OUTPUT_FILE}
