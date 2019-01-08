# FROM ensemblorg/ensembl-vep:release_93
FROM ensemblorg/ensembl-vep@sha256:abc7a41b218c4e051bf97ed5e42add0b3c9fcf84929da28f1f80a215f259a69e

USER root

RUN apt-get update && apt-get install -y \
    git \
    python-pip \
    python-dev \
    libmysqlclient-dev

WORKDIR /app

# pyhgvs 0.9.4
RUN pip install git+git://github.com/counsyl/hgvs.git@aebe5bd9683f4b5937fd653ce4b13fcd4f3ebb10
# maxentpy 0.0.1
RUN pip install git+https://github.com/kepbod/maxentpy.git

# RUN pip install --upgrade pip NOTE: pyhgvs breaks with pip > v8
ADD requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

ADD . /app

ENTRYPOINT ["python", "calcVarPriors.py"]
