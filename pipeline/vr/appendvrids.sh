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


# ...and wait for the HTTP endpoint to be ready (not just TCP port open)
MAX_WAIT=300
i=0
until curl -sf http://localhost:5000/seqrepo/ping > /dev/null 2>&1; do
    i=$((i + 1))
    if [ $i -ge $MAX_WAIT ]; then
        echo "Timed out waiting for seqrepo-rest-service to become ready after ${MAX_WAIT}s" >&2
        exit 1
    fi
    # Fail fast if the container has already exited
    if ! docker ps --filter name=seqrepo-rest-service --filter status=running --format '{{.Names}}' | grep -q seqrepo-rest-service; then
        echo "seqrepo-rest-service container exited unexpectedly (check: docker logs seqrepo-rest-service)" >&2
        exit 1
    fi
    sleep 1
done
echo " => ...seqrepo-rest-service is ready"


# ----------------------------------------------------
# --- 2. execute append-vr-ids
# ----------------------------------------------------
PATH=../utilities:${PATH} python3 appendVRIds.py \
    -i ${ARTIFACT_DIR}/${INPUT_FILE} \
    -o ${ARTIFACT_DIR}/${OUTPUT_FILE}
