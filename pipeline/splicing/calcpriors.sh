#!/usr/bin/env bash
set -o nounset
set -o errexit

REFERENCES_DIR=$1
INPUT_OUTPUT_DIR=$2
INPUT_FILE=$3
OUTPUT_FILE=$4

# Download and verify references - to be called from within the docker
echo "Downloading references. If not already present, this may take 1-2 hours..."
docker run -it --rm \
    --user=`id -u`:`id -g` \
    -v ${REFERENCES_DIR}:/references \
    brcachallenge/splicing-pipeline references
echo "Reference download and installation complete."

# Run short test to ensure proper setup
echo "Running short test to ensure proper set up."
# TODO: exit if tests error
docker run -it --rm \
    --user=`id -u`:`id -g` \
    -v ${REFERENCES_DIR}:/references:ro \
    brcachallenge/splicing-pipeline test short
echo "Tests passed!"

# Calculate priors
echo "Calculating priors, this may take a few hours..."
docker run --rm -it \
    --user=`id -u`:`id -g` \
    -v ${REFERENCES_DIR}:/references:ro \
    -v ${INPUT_OUTPUT_DIR}:/data \
    brcachallenge/splicing-pipeline calc ${INPUT_FILE} ${OUTPUT_FILE}
echo "Prior calculation complete!"

