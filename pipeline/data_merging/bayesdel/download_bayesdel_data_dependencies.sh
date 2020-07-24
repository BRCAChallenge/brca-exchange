#!/usr/bin/env bash
set -e

# downloads the data dependency for the bayesdel part of the pipeline

if [[ "$#" -ne 1 ]]; then
    echo "Usage: $0 [data_dir]"
    exit 1
fi

DATA_DIR=$1
DONE_FILE="${DATA_DIR}/DONE"

echo "Working in ${DATA_DIR}"
mkdir -p ${DATA_DIR}

if [[ -f "${DONE_FILE}" ]]; then
    echo "Data seems already downloaded in ${DATA_DIR}. Skipping download."
    exit 0
fi

cd ${DATA_DIR}

echo "Downloading download scripts..."
# saving the download scripts along with the data, such that the provenance of the data can be better tracked.
curl -L http://www.bjfenglab.org/download/get_GRCh38_inst.sh  > get_GRCh38_inst.sh
curl -L http://www.bjfenglab.org/download/get_GRCh38_update.sh > get_GRCh38_update.sh

echo "Downloading main data..."
bash get_GRCh38_inst.sh
echo "Downloading data updates..."
bash get_GRCh38_update.sh

date > ${DONE_FILE}
