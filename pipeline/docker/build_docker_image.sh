#!/bin/bash

# since we need this file during image building, we need to stage it into the docker context
cp ../requirements.txt requirements_docker.txt

docker build -t brca-exchange-pipeline .

rm requirements_docker.txt

