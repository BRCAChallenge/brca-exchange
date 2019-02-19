#!/bin/bash

REF_PATH="/Users/Faisal/Resources/brca-exch/resources"
DATA_PATH="/Users/Faisal/Resources/brca-exch/past-releases/release-11-03-18/output/release/artifacts"

# --user=`id -u`:`id -g` \

#docker run --rm -it \
#	--entrypoint /bin/bash \
#	-v "${REF_PATH}":/references \
#	-v "${DATA_PATH}":/data \
#	-v `pwd`:/app \
#	brcachallenge/splicing-pipeline

LAST_COMMIT=$( git log -n 1 --pretty=format:%h -- . )
DOCKER_IMG_ID="brcachallenge/splicing-pipeline/${LAST_COMMIT}"

docker build -t ${DOCKER_IMG_ID}

docker run -it --rm --user=`id -u`:`id -g` \
  -v "${REF_PATH}":/references:ro \
  -v "${DATA_PATH}":/data \
  -v "$(pwd)":/app \
  ${DOCKER_IMG_ID} \
  calc /data/built_with_mupit.tsv /data/built_with_priors_compare.tsv
