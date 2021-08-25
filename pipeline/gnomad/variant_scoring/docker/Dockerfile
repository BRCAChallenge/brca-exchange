FROM continuumio/miniconda3

RUN mkdir /opt/variant_scoring
WORKDIR /opt/variant_scoring

COPY . .

RUN apt-get update && apt-get install -y build-essential
RUN conda env update --name base --file gnomad/variant_scoring/environment.yml


