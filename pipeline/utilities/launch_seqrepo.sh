#!/usr/bin/env bash

SEQ_REPO_DIR="${1:-/usr/local/share/seqrepo/latest}"

#
# Launch the seqrepo rest API docker container
docker rm -f seqrepo-rest-service 2>/dev/null || true
docker run --name seqrepo-rest-service \
       --detach -p 5000:5000 \
       --memory 8g \
       --restart on-failure \
       -v ${SEQ_REPO_DIR}:/mnt/seqrepo \
       biocommons/seqrepo-rest-service \
       seqrepo-rest-service /mnt/seqrepo
