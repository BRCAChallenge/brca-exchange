#!/usr/bin/env bash

set -e

if [[ "$#" -lt 2 ]]; then
    echo "$#"
    echo "Usage: $0 vcf_file victor_data_dir wdir [image]"
    exit 1
fi

INPUT_VCF_HOST=$1
VICTOR_DATA_DIR_HOST=$2
WDIR_HOST=$3
IMAGE="${4:-brcachallenge/victor:0.1}"

INPUT_FILE_DOCKER=/tmp/my_input.vcf
WDIR_DOCKER=/opt/victor/wdir

echo "running in image ${IMAGE}"

docker run --rm \
  --user=`id -u`:`id -g` \
  -v ${VICTOR_DATA_DIR_HOST}:/opt/victor/VICTOR/data:ro \
  -v ${INPUT_VCF_HOST}:${INPUT_FILE_DOCKER}:ro \
  -v ${WDIR_HOST}:${WDIR_DOCKER} \
  ${IMAGE} ${INPUT_FILE_DOCKER} ${WDIR_DOCKER}
