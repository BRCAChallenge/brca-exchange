# FROM ubuntu:16.04
FROM ensemblorg/ensembl-vep:release_93

USER root

# RUN perl /opt/vep/src/ensembl-vep/INSTALL.pl

RUN apt-get update && apt-get install -y \
    git \
    python-pip \
    python-dev \
    libmysqlclient-dev

WORKDIR /app

RUN pip install --upgrade pip
ADD requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# pyhgvs 0.9.4
RUN pip install git+git://github.com/counsyl/hgvs.git@aebe5bd9683f4b5937fd653ce4b13fcd4f3ebb10

# USER vep

ADD . /app

ENTRYPOINT ["python", "calcVarPriors.py", "-b", "enigma", "-v", "mod_res_dn_brca20160525.txt"]
