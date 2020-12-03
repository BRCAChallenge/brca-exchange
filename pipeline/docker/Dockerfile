FROM python:3.8

RUN chmod 1777 /tmp

RUN apt-get update && apt-get install -y \
    curl \
    git \
    libmariadb-dev-compat \
    libmariadb-dev \
    liblzo2-dev \
    libbz2-dev \
    liblzma-dev \
    pkg-config \
    python3 \
    python3-pip \
    python3-gdbm \
    python-lzo \
    rsync \
    vim \
    wget \
    zlib1g-dev \
    netcat

# Get the Docker binary
RUN curl -fsSL get.docker.com -o get-docker.sh \
    && sh get-docker.sh

WORKDIR /opt

COPY pipeline/requirements.txt .
COPY test-requirements.txt .

# pip 20.3+ uses strict dependency resolver that causes biocommons/bioutils and hgvs/ipython errors
RUN pip install pip==20.2

# install numpy first to avoid issues with bio python and bx-python (see also https://github.com/LUMC/vep2lovd/issues/1)
RUN pip install $(grep numpy requirements.txt)

RUN pip install $(grep -i cython requirements.txt)
RUN pip install -r requirements.txt -r test-requirements.txt

# install vcf tools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
RUN tar zxf vcftools*.tar.gz
RUN cd vcftools* && ./configure && make && make install && cd .. && rm -r vcftools*

# install tabix
RUN wget https://downloads.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2
RUN tar jxf tabix*.tar.bz2
RUN cd tabix* && make && cp tabix /usr/local/bin && cd .. && rm -r tabix*

ARG res=/files/resources

ENV BRCA_RESOURCES=$res

RUN mkdir -p $res /files/data && chmod -R o+rwx /files

RUN rm -r  /root/.cache

ARG FORCE_REBUILD=0
COPY . /opt/brca-exchange

ENV LUIGI_CONFIG_PATH="/opt/luigi_pipeline_credentials.cfg"

ARG IS_GIT_DIRTY="False"
ARG GIT_COMMIT=""
LABEL GitCommit=$GIT_COMMIT IsGitDirty=$IS_GIT_DIRTY

CMD ["/opt/brca-exchange/pipeline/docker/run_luigi.sh"]

