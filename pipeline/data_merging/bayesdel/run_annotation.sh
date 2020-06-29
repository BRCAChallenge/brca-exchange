#!/usr/bin/env bash
set -e

if [[ "$#" -ne 2 ]]; then
    echo "got argumetns $#"
    echo "Usage: $0 [INPUT_FILE] [WORK_DIR]"
    exit 1
fi

echo "calling run annotation $1 and $2"
ORIG_FILE=$1
WDIR=$2

VCF_FILE=${WDIR}/input.vcf

cd ${WDIR}

echo "copying ${ORIG_FILE}"
echo $(ls -l ${ORIG_FILE})

# prepare files
cp ${ORIG_FILE} ${VCF_FILE}
bgzip -f ${VCF_FILE}
tabix -f ${VCF_FILE}.gz

# run annotation
export SLURM_ARRAY_TASK_ID=1
export VICTOR_GENOME=GRCh38

export VCF=${VCF_FILE}.gz
export OUT=output

# TODO: limit threads?
slurm.annotate
