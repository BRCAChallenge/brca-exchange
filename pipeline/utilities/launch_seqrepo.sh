#!/usr/bin/env bash

SEQ_REPO_DIR="${1:-/usr/local/share/seqrepo/latest}"

#
# Launch the seqrepo rest API docker container
docker run --name seqrepo-rest-service \
       --detach --rm -p 5000:5000 \
       -v ${SEQ_REPO_DIR}:/mnt/seqrepo \
       biocommons/seqrepo-rest-service \
       /mnt/seqrepo
