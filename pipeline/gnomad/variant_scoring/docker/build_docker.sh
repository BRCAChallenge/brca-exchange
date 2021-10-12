#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd ../../..
pwd
docker build --file gnomad/variant_scoring/docker/Dockerfile  -t variant_scoring .

