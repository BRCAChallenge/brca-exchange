#!/usr/bin/env bash

# (paths and configs with respect to dev machine)
DATA_DIR=/data/variant_scoring/data

# the DATA_DIR is expected to contain the following
#├── original_gnomad_data                                      # downloaded e.g. using download_data.sh
#│   ├── gnomad.exomes.coverage.summary.tsv.bgz
#│   ├── gnomad.exomes.r2.1.1.sites.13.vcf.bgz
#│   ├── gnomad.exomes.r2.1.1.sites.17.vcf.bgz
#│   ├── gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz
#│   ├── gnomad.genomes.v3.1.1.sites.chr13.vcf.bgz
#│   └── gnomad.genomes.v3.1.1.sites.chr17.vcf.bgz
#└── resources                                                 # taken e.g. from main brca-exchange pipeline
#    ├── hg19ToHg38.over.chain.gz
#    ├── hg38.fa


UTA_PORT=50828
UTA_RELEASE_DATE=20171026
SEQ_REPO_DIR_DOCKER=/mnt/seq_repo

COMMON_DOCKER_OPT=" --rm -v ${DATA_DIR}:/mnt --network host -v /data/variant_scoring/code/pipeline:/opt/variant_scoring -v /data/seqrepo:${SEQ_REPO_DIR_DOCKER} -e UTA_DB_URL=postgresql://anonymous@0.0.0.0:${UTA_PORT}/uta/uta_${UTA_RELEASE_DATE} -e HGVS_SEQREPO_DIR=${SEQ_REPO_DIR_DOCKER}/2021-01-29"

# run process_gnomad_data.py
docker run ${COMMON_DOCKER_OPT} variant_scoring bash -c "python /opt/variant_scoring/gnomad/variant_scoring/process_gnomad_data.py --cores 16 --mem-per-core-mb 2400 --gene-config-path /opt/variant_scoring/workflow/gene_config_brca_only.txt /mnt/original_gnomad_data /mnt/processed_brca"

# run varint_scoring.py
docker run ${COMMON_DOCKER_OPT} variant_scoring bash -c "python /opt/variant_scoring/gnomad/variant_scoring/variant_scoring.py --gene-config-path /opt/variant_scoring/workflow/gene_config_brca_only.txt --resource-dir /mnt/resources /mnt/processed_brca /mnt/variant_scoring_df.parquet"

