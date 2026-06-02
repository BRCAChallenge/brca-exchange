#!/usr/bin/env bash

DOCKER_IMAGE_NAME="brcachallenge/exlovd-extraction:0.1"

INPUT_OUTPUT_DIR=$1
echo "dir: ${INPUT_OUTPUT_DIR}"
shift

echo "running with $@"

docker run --rm \
    --user=`id -u`:`id -g` \
    -v ${INPUT_OUTPUT_DIR}:/data \
    ${DOCKER_IMAGE_NAME} $@
