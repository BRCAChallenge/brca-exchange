#!/bin/bash

# note: all paths are valid within the docker container!

# if DATA_DATE is set as environment variable, set the pipeline date explicitely (useful for continuing a pipeline run on another day)
if [ -n "${DATA_DATE}" ]; then
    DATE_PARAM_OPT="--date ${DATA_DATE}"
fi

OUTPUT_DIR=/files/data/output
PARENT_DIR=/files/data
BRCA_RESOURCES=/files/resources

PREVIOUS_RELEASE_TAR=/files/previous_release.tar.gz

RELEASE_NOTES=/files/release_notes.txt

CODE_MNT=$(mount | grep /opt/brca-exchange)
[ -z "${CODE_MNT}" ] || echo "WARNING: BRCA Code base mounted from host file system"

cd /opt/brca-exchange

echo "Running brca exchange pipeline:"
echo "Git hash: $(git log | head -n 1)"

cd /opt/brca-exchange/pipeline/luigi

python -m luigi --logging-conf-file luigi_log_configuration.conf --module CompileVCFFiles RunAll --resources-dir ${BRCA_RESOURCES} --file-parent-dir ${PARENT_DIR} --output-dir ${OUTPUT_DIR} --previous-release-tar ${PREVIOUS_RELEASE_TAR} --release-notes ${RELEASE_NOTES} ${DATE_PARAM_OPT} --local-scheduler
