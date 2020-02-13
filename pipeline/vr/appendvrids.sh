#!/bin/bash

ARTIFACT_DIR="$1"
INPUT_FILE="$2"
OUTPUT_FILE="$3"
IMAGE_NAME="${4:-brcachallenge/append-vr-ids:0.1}"
SEQ_REPO_DIR="${5:-/usr/local/share/seqrepo}"

# --- pre-step: build vr image (caching will make this fast if it was already built)
docker build -t "${IMAGE_NAME}" .

# ----------------------------------------------------
# --- 1. spin up seqrepo-rest-service
# ----------------------------------------------------
docker run --rm -d \
  --user=`id -u`:`id -g` \
  --name seqrepo-service \
  --network=host \
  -v "${SEQ_REPO_DIR}":/usr/local/share/seqrepo \
  biocommons/seqrepo-rest-service:latest

# ...and wait for it to be available
./wait-for localhost:5000 && echo " => ...seqrepo-rest-service is ready"

# ----------------------------------------------------
# --- 2. execute append-vr-ids
# ----------------------------------------------------
docker run --rm \
  --name append-vr-ids \
  --network=host \
  -v "${ARTIFACT_DIR}":/artifacts \
  -e UTA_DB_URL \
  "${IMAGE_NAME}" \
  -i "/artifacts/${INPUT_FILE}" \
  -o "/artifacts/${OUTPUT_FILE}"

# --- post-step: when vr's done, bring down its dependent services
docker stop seqrepo-service
