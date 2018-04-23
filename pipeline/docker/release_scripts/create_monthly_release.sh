#!/bin/bash

# Script to create data for a monthly release. Sets up environment and kicks off pipeline run

if [ "$#" -lt "2" ] || [ "$#" -gt "3" ]; then
    echo "Arguments:"
    echo " 1. git commit hash -- hash of the main git repository from which to use the code"
    echo " 2. path to previous release tar file -- used to calculate the delta between the release"
    echo " 3. (optional) release tag in the format of data_release_yyyy-MM-dd"
    exit 1
fi

function msg
{
    echo "$(date) $1"
}

DATA_DATE="$(date +%Y-%m-%d)"

GIT_COMMIT=$1
PREVIOUS_RELEASE=$2

DEFAULT_RELEASE_TAG="data_release_${DATA_DATE}"
RELEASE_TAG=${3:-${DEFAULT_RELEASE_TAG}} # for docker image and tag in git

DEFAULT_GIT_REPO="https://github.com/BRCAChallenge/brca-exchange.git"
BRCA_GIT_REPO=${BRCA_GIT_REPO:-${DEFAULT_GIT_REPO}}

BASE_DIR="/home/pipeline/monthly_releases"
WORK_DIR="/home/pipeline/monthly_releases/${RELEASE_TAG}"
LOG_DIR="${BASE_DIR}/logs"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -d "${WORK_DIR}" ]; then
    msg "Directory ${WORK_DIR} already exists!"
    exit 1
fi

RELEASE_NOTES="${BASE_DIR}/release_notes/release_notes_${RELEASE_TAG}"

if [ ! -f "${RELEASE_NOTES}" ]; then
    msg "WARNING: files ${RELEASE_NOTES} doesn't exist. Creating it" 
    touch ${RELEASE_NOTES}
fi

mkdir -p ${WORK_DIR}/brca_out

# Cloning current repo
msg "Setting up git repository"

CODE_BASE=${WORK_DIR}/brca-exchange
cd ${WORK_DIR}

git clone ${BRCA_GIT_REPO}
cd ${CODE_BASE}

git checkout ${GIT_COMMIT}

# Create docker image
msg "Creating docker image"

cd ${CODE_BASE}/pipeline/docker
./build_docker_image.sh ${GIT_COMMIT} ${RELEASE_TAG} || { msg "Something wrong with docker image creation"; exit 1; }

# Download reference data
msg "Download reference data..."
RESOURCES=${WORK_DIR}/resources
mkdir -p ${RESOURCES}

${CODE_BASE}/pipeline/download_resources_files.sh ${RESOURCES} >> "${LOG_DIR}/download_resources_${RELEASE_TAG}.log" 2>&1 || { msg "Downloading resource files failed"; exit 1; }

# run pipeline
msg "Kicking off pipeline!"

PIPELINE_BASE_COMMAND="${BASE_DIR}/scripts/run_docker.sh ${WORK_DIR} brcachallenge/brca-exchange-pipeline:${RELEASE_TAG} ${PREVIOUS_RELEASE} ${RELEASE_NOTES} ${DATA_DATE}"

PIPELINE_COMMAND="${PIPELINE_BASE_COMMAND} > ${LOG_DIR}/pipeline_run_${RELEASE_TAG}_\$(date +%Y%m%d_%H%M%S).log 2>&1"

POSTPROCESSING_COMMANDS_FILE="${SCRIPT_DIR}/postprocessing_commands/postprocessing_commands_${RELEASE_TAG}.sh"

# requires pip install j2cli
PIPELINE_COMMAND=${PIPELINE_COMMAND} RELEASE_NOTES=${RELEASE_NOTES} RELEASE_TAG=${RELEASE_TAG} j2 ${SCRIPT_DIR}/postprocessing_commands_template.j2 > "${POSTPROCESSING_COMMANDS_FILE}"


msg "In case of errors, restart the pipeline using the command below"
msg "(See also ${POSTPROCESSING_COMMANDS_FILE})"
echo ""
echo "${PIPELINE_COMMAND}"

eval "${PIPELINE_COMMAND}"

printf ' \n\n'
msg "Luigi says:"
grep 'This progress looks' ${LOG_DIR}/pipeline_run_${RELEASE_TAG}*.log
printf ' \n\n'

echo "Running variantsBySource.py script:"
${PIPELINE_BASE_COMMAND} 'python /opt/brca-exchange/pipeline/utilities/variantsBySource.py  -i /files/data/output/release/built_with_change_types.tsv -c true'

msg "If data seems okay, you can find the next steps in ${POSTPROCESSING_COMMANDS_FILE}"
