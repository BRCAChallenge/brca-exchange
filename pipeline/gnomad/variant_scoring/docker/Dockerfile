FROM continuumio/miniconda3

RUN mkdir -p /opt/variant_scoring /data/resources
WORKDIR /opt/variant_scoring


COPY . .

# workaround to install Java runtime
RUN mkdir -p /usr/share/man/man1/ && touch /usr/share/man/man1/java.1.gz.dpkg-tmp
RUN apt-get update && apt-get install -y build-essential default-jre libcurl4-openssl-dev zlib1g-dev
RUN conda env update --name base --file gnomad/variant_scoring/environment.yml

ENV PYTHONPATH=/opt/variant_scoring


