FROM ensemblorg/ensembl-vep:release_93

USER root

RUN apt-get update && apt-get install -y \
    git \
    python-pip \
    python-dev \
    libmysqlclient-dev

WORKDIR /app

# pyhgvs 0.9.4 - NOTE: Will not work with newer pip
RUN pip install git+git://github.com/counsyl/hgvs.git@aebe5bd9683f4b5937fd653ce4b13fcd4f3ebb10

# RUN pip install --upgrade pip
ADD requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

ADD . /app

ENTRYPOINT ["python", "calcVarPriors.py"]
