#!/bin/bash

# Creates a docker image to run the BRCA Exchange data pipeline
#
# Optional arguments:
#  1. tag for the docker image
#
# The codebase is included as-is into the docker image. Note, that in this process all files not excluded in the '.dockerignore' file get copied into the image, which may significantly slow down the building process.
# If you want to build a docker image with a specific branch or commit, check it out before calling this script.

DEFAULT_DOCKER_TAG="latest"
DOCKER_TAG=$1

DOCKER_TAG=${DOCKER_TAG:-${DEFAULT_DOCKER_TAG}}

COMMIT_LENGTH=${#COMMIT}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${SCRIPT_DIR}

# checking whether we are in dirty git state (ignoring untracked files though)
IS_GIT_DIRTY="False"
if [  "$(git status --short | grep -cv '^??')" -gt 0 ]; then
    echo "WARNING: files have been modified! Please consider committing your changes."
    IS_GIT_DIRTY="True"
fi

BUILD_DIR="$(realpath ${SCRIPT_DIR}/../..)"

echo "Build context is ${BUILD_DIR}"

DOCKER_IMAGE_NAME="brcachallenge/brca-exchange-pipeline:${DOCKER_TAG}"
COMMIT=$(git rev-parse HEAD)

echo "Building ${DOCKER_IMAGE_NAME} with ${COMMIT}"
docker build -f "${BUILD_DIR}/pipeline/docker/Dockerfile" -t ${DOCKER_IMAGE_NAME} --build-arg FORCE_REBUILD=1 --build-arg IS_GIT_DIRTY=${IS_GIT_DIRTY} --build-arg GIT_COMMIT=${COMMIT} ${BUILD_DIR}

DOCKER_EXIT="$?"

exit ${DOCKER_EXIT}
