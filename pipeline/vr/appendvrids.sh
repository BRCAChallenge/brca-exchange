#!/bin/bash

INPUT_DIR="$1"
OUTPUT_FILE="$2"
IMAGE_NAME="${3:-brcachallenge/append-vr-ids:0.1}"
SEQ_REPO_DIR="${4:-/usr/local/share/seqrepo}"

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
  -v "${INPUT_DIR}":/output \
  -e UTA_DB_URL \
  "${IMAGE_NAME}" \
  -i "/output" \
  -a "/output/release/artifacts" \
  -o "/output/release/artifacts/${OUTPUT_FILE}"
