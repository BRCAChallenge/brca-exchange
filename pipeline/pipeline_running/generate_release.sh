#!/bin/bash

# TODO: add command line wrapper
# TODO: crontab: git pull master

ROOT_DIR=$(realpath $1)
CREDENTIALS_PATH=$(realpath $2)
PREVIOUS_RELEASE_DIR=$(realpath $3)

GIT_COMMIT=${4-:"master"}

DATA_DATE="$(date +%Y-%m-%d)"
WDIR="${ROOT_DIR}/data_release_${DATA_DATE}"
mkdir -p ${WDIR}


CODE_BASE="${WDIR}/code"
[ -d ${CODE_BASE} ] || git clone https://github.com/BRCAChallenge/brca-exchange.git ${CODE_BASE}
cd ${CODE_BASE} && git checkout ${GIT_COMMIT}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

CFG_PATH=${CODE_BASE}/pipeline/brca_pipeline_cfg.mk

jinja2 -D DATA_DATE="${DATA_DATE}" -D WORK_DIR="${WDIR}" \
       -D CODE_BASE="${CODE_BASE}" \
       -D CREDENTIALS_PATH="${CREDENTIALS_PATH}" \
       -D PREVIOUS_RELEASE_DIR="${PREVIOUS_RELEASE_DIR}" \
       -D GIT_COMMIT="${GIT_COMMIT}" \
       "${CODE_BASE}/pipeline/pipeline_running/brca_pipeline_cfg.mk.j2" > ${CFG_PATH}

cd "${CODE_BASE}/pipeline"

make build-release

echo "You can issue pipeline commands using"
echo "make CONFIG_PATH=${CFG_PATH} [cmd]"
echo " -- or"
echo "cd ${CODE_BASE}/pipeline && make [cmd]"
