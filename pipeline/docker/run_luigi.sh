#!/bin/bash

# note: all paths are valid within the docker container!

OUTPUT_DIR=/files/data/output
PARENT_DIR=/files/data
BRCA_RESOURCES=/files/resources

PREVIOUS_RELEASE=/files/previous_release/built_with_change_types.tsv
RELEASE_NOTES=/files/release_notes.txt

cd /opt/brca-exchange

echo "Running brca exchange pipeline:"
echo "Git hash: $(git log | head -n 1)"

cd /opt/brca-exchange/pipeline/luigi

python -m luigi --logging-conf-file luigi_log_configuration.conf --module CompileVCFFiles RunAll --resources-dir ${BRCA_RESOURCES} --file-parent-dir ${PARENT_DIR} --output-dir ${OUTPUT_DIR} --previous-release ${PREVIOUS_RELEASE} --previous-release-date ${PREVIOUS_RELEASE_DATE} --release-notes ${RELEASE_NOTES} --local-scheduler
