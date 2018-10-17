#!/usr/bin/env bash
set -o nounset
set -o errexit

REFERENCES_DIR=$1
INPUT_OUTPUT_DIR=$2
INPUT_FILE_NAME=$3
OUTPUT_FILE_NAME=$4
DOCKER_IMAGE_NAME=${5:-brcachallenge/splicing-pipeline}

# Download and verify references - to be called from within the docker
echo "Downloading references. If not already present, this may take 1-2 hours..."
docker run --rm \
    --user=`id -u`:`id -g` \
    -v ${REFERENCES_DIR}:/references \
    ${DOCKER_IMAGE_NAME} references
echo "Reference download and installation complete."

# Run short test to ensure proper setup
echo "Running short test to ensure proper set up."
docker run --rm \
    --user=`id -u`:`id -g` \
    -v ${REFERENCES_DIR}:/references:ro \
    ${DOCKER_IMAGE_NAME} test short
echo "Tests passed!"

# Calculate priors
echo "Calculating priors, this may take a few hours..."
docker run --rm \
    --user=`id -u`:`id -g` \
    -v ${REFERENCES_DIR}:/references:ro \
    -v ${INPUT_OUTPUT_DIR}:/data \
    ${DOCKER_IMAGE_NAME} --processes $(nproc) calc "/data/${INPUT_FILE_NAME}" "/data/${OUTPUT_FILE_NAME}"
echo "Prior calculation complete!"

