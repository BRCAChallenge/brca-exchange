#!/usr/bin/env bash
set -o nounset
set -o errexit

INPUT_OUTPUT_DIR=$1
INPUT_FILE_NAME=$2
OUTPUT_FILE_NAME=$3

DOCKER_IMAGE_NAME=${5:-brcachallenge/append-vr-ids}

# Get VR Id's
echo "Gathering VR Id's..."
docker run --rm \
    --user=`id -u`:`id -g` \
    -v ${INPUT_OUTPUT_DIR}:/data \
    ${DOCKER_IMAGE_NAME} \
    -i /data/${INPUT_FILE_NAME} \
    -o /data/${OUTPUT_FILE_NAME}
echo "Done!"

