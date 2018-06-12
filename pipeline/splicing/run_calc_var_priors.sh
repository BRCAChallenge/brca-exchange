#!/bin/bash

# note: all paths are valid within the docker container!

# if DATA_DATE is set as environment variable, set the priors date explicitely (useful for continuing a priors run on another day)
if [ -n "${DATA_DATE}" ]; then
    DATE_PARAM_OPT="--date ${DATA_DATE}"
fi

REFSEQ_ANNOTATION=/files/resources/refseq_annotation.hg38.gp
HG38=/files/resource/hg38.fa
INPUT=/files/data/output/release/artifacts/built.tsv
OUTPUT=/files/data/output/release/artifacts/built_with_priors.tsv


CODE_MNT=$(mount | grep /opt/brca-exchange)
[ -z "${CODE_MNT}" ] || echo "WARNING: BRCA Code base mounted from host file system"

cd /opt/brca-exchange

echo "Running brca exchange pipeline:"
echo "Git hash: $(git log | head -n 1)"

cd /opt/brca-exchange/pipeline/splicing

python calcVarPriors.py -i ${INPUT} -o ${OUTPUT} -b "enigma" -v mod_res_dn_brca20160525.txt -g ${HG38} -t ${REFSEQ_ANNOTATION}