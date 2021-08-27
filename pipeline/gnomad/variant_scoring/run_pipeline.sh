#!/usr/bin/env bash

# (paths and configs with respect to dev machine)
DATA_DIR=/data/variant_scoring/data

COMMON_DOCKER_OPT=" --rm -v ${DATA_DIR}:/mnt -v /data/variant_scoring/code/pipeline:/opt/variant_scoring" 

docker run ${COMMON_DOCKER_OPT} variant_scoring bash -c "python /opt/variant_scoring/gnomad/variant_scoring/process_gnomad_data.py --cores 16 --mem-per-core-mb 2400 --gene-config-path /opt/variant_scoring/workflow/gene_config_brca_only.txt /mnt/original_gnomad_data /mnt/processed_brca"

docker run ${COMMON_DOCKER_OPT} variant_scoring bash -c "python /opt/variant_scoring/gnomad/variant_scoring/variant_scoring.py  --resource-dir /mnt/resources /mnt/processed_brca /mnt/variant_scoring_df.parquet"

