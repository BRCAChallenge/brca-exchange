#!/usr/bin/env bash

SEQ_REPO_DIR="${1:-/usr/local/share/seqrepo/latest}"

# If given a parent seqrepo directory rather than an instance, resolve to latest
if [ ! -f "${SEQ_REPO_DIR}/aliases.sqlite3" ] && [ -d "${SEQ_REPO_DIR}/latest" ]; then
    SEQ_REPO_DIR="${SEQ_REPO_DIR}/latest"
fi

#
# Launch the seqrepo rest API docker container
docker run --name seqrepo-rest-service \
       --detach --rm -p 5000:5000 \
       -v ${SEQ_REPO_DIR}:/mnt/seqrepo \
       biocommons/seqrepo-rest-service \
       seqrepo-rest-service /mnt/seqrepo
