#!/bin/bash

# Creates a docker image to run the BRCA Exchange data pipeline
#
# Optional arguments:
#  1. git commit: either a branch name or a hash
#  2. tag for the docker image
#
# If no arguments are provided, the code of the current master branch is used and the docker imaged tagged with "latest"

DEFAULT_GIT_REPO="https://github.com/BRCAChallenge/brca-exchange.git"
DEFAULT_COMMIT="master"
DEFAULT_DOCKER_TAG="latest"

SHA1_HASH_HEX_LENGTH=40

COMMIT=$1
DOCKER_TAG=$2

BRCA_GIT_REPO=${BRCA_GIT_REPO:-${DEFAULT_GIT_REPO}}
COMMIT=${COMMIT:-${DEFAULT_COMMIT}}
DOCKER_TAG=${DOCKER_TAG:-${DEFAULT_DOCKER_TAG}}

COMMIT_LENGTH=${#COMMIT}

if [ "${COMMIT_LENGTH}" -eq "${SHA1_HASH_HEX_LENGTH}" ]; then
    FORCE_REBUILD=0
else
    echo "No commit hash given. Will force to clone brca-exchange repository and checkout branch ${COMMIT}"
    FORCE_REBUILD=$(date +%s)
fi

# since we need this file during image building, we need to stage it into the docker context
cp ../requirements.txt requirements_docker.txt

DOCKER_IMAGE_NAME="brcachallenge/brca-exchange-pipeline:${DOCKER_TAG}"
echo "Building ${DOCKER_IMAGE_NAME} with code from ${BRCA_GIT_REPO} ${COMMIT}"
docker build -t ${DOCKER_IMAGE_NAME} --build-arg FORCE_REBUILD=${FORCE_REBUILD} --build-arg BRCA_GIT_REPO=${BRCA_GIT_REPO} --build-arg BRCA_EXCHANGE_COMMIT=${COMMIT} .

DOCKER_EXIT="$?"

rm requirements_docker.txt

exit ${DOCKER_EXIT}
