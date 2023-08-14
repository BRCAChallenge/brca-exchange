#!/bin/bash

ARTIFACT_DIR="$1"
INPUT_FILE="$2"
OUTPUT_FILE="$3"
SEQ_REPO_DIR="${4:-/usr/local/share/seqrepo}"

echo "Hello World"
echo ${SEQ_REPO_DIR}

# --- pre-step: build the seqrepo-rest-service docker
# See https://github.com/biocommons/seqrepo-rest-service

# ----------------------------------------------------
# --- 1. spin up seqrepo-rest-service
# ----------------------------------------------------
docker run \
  --name seqrepo-rest-service \
  --user=`id -u`:`id -g` \
  --detach --rm -p 5000:5000 \
  --network=host \
  -v ${SEQ_REPO_DIR}:/usr/local/share/seqrepo \
  biocommons/seqrepo-rest-service \
  seqrepo-rest-service /usr/local/share/seqrepo


# ...and wait for it to be available
./wait-for localhost:5000 && echo " => ...seqrepo-rest-service is ready"


# ----------------------------------------------------
# --- 2. execute append-vr-ids
# ----------------------------------------------------
python3 appendVRIds.py -i ${ARTIFACT_DIR}/${INPUT_FILE} \
       -o ${ARTIFACT_DIR}/${OUTPUT_FILE}
