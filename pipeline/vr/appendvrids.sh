#!/bin/bash

ARTIFACT_DIR="$1"
INPUT_FILE="$2"
OUTPUT_FILE="$3"
IMAGE_NAME="$4"

# pre-step: build vr image (caching will make this fast if it was already built)
docker build -t "${IMAGE_NAME}" .

# spin up seqrepo-rest-service...
docker run --rm -d \
  --name seqrepo-service \
  -p 5000:5000 --network=host \
  -v /files/resources/seq_repo:/usr/local/share/seqrepo \
  biocommons/seqrepo-rest-service:latest
# ...and wait for it to be available
./wait-for localhost:5000

# now we can run append-vr-ids...
docker run --rm \
  --name append-vr-ids \
  --network=host \
  -v "$ARTIFACT_DIR":/artifacts \
  -e UTA_DB_URL \
  "${IMAGE_NAME}" \
  -i "/artifacts/${INPUT_FILE}" \
  -o "/artifacts/${OUTPUT_FILE}"
# ...and wait for *it* to complete
docker wait append-vr-ids

# and when vr's done, bring down its dependent services
docker stop seqrepo-service
