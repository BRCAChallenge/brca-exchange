#!/bin/bash

# if DATA_DATE is set as environment variable, set the pipeline date explicitely (useful for continuing a pipeline run on another day)
if [ -n "${DATA_DATE}" ]; then
    DATE_PARAM_OPT="--date ${DATA_DATE}"
fi

# note: paths are valid within the docker container!
PARENT_DIR=/files/data
OUTPUT_DIR="${PARENT_DIR}/output"
BRCA_RESOURCES=/files/resources

# TODO doc why different mechanism
if [ "$#" -ne "3" ]; then
    echo "Require host (!) paths for references directory and data directory"
    echo "Usage: run_luigi.sh [PRIORS_REFERENCES] [OUTPUT_DIR_HOST] [PRIORS_DOCKER_IMAGE_NAME] [TASK]" 
fi

PRIORS_REFERENCES=$1
OUTPUT_DIR_HOST="$2/output"
PRIORS_DOCKER_IMAGE=$3
LUIGI_TASK=${4:-"RunAll"}

PREVIOUS_RELEASE_TAR=/files/previous_release.tar.gz

RELEASE_NOTES=/files/release_notes.txt

CODE_MNT=$(mount | grep /opt/brca-exchange)
[ -z "${CODE_MNT}" ] || echo "WARNING: BRCA Code base mounted from host file system"

cd /opt/brca-exchange

echo "Running brca exchange pipeline:"
echo "Git hash: $(git log | head -n 1)"

cd /opt/brca-exchange/pipeline/luigi

echo "Attempting to run task ${LUIGI_TASK}"
python -m luigi --logging-conf-file luigi_log_configuration.conf --module CompileVCFFiles ${LUIGI_TASK} --resources-dir ${BRCA_RESOURCES} --file-parent-dir ${PARENT_DIR} --output-dir ${OUTPUT_DIR} --previous-release-tar ${PREVIOUS_RELEASE_TAR} --priors-references-dir ${PRIORS_REFERENCES} --priors-docker-image-name ${PRIORS_DOCKER_IMAGE} --output-dir-host ${OUTPUT_DIR_HOST} --release-notes ${RELEASE_NOTES} ${DATE_PARAM_OPT} --local-scheduler
