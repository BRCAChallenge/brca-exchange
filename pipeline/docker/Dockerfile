FROM python:3.12-bookworm

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
    python3-gdbm \
    rsync \
    vim \
    wget \
    zlib1g-dev \
    netcat-openbsd

# Get the Docker binary
RUN curl -fsSL get.docker.com -o get-docker.sh \
    && sh get-docker.sh

WORKDIR /opt

COPY pipeline/requirements.txt .
COPY test-requirements.txt .

# pip 20.3+ uses strict dependency resolver that causes biocommons/bioutils and hgvs/ipython errors
# Using pip 24.3.1 for Python 3.12 compatibility
RUN pip install --upgrade pip==24.3.1

# install six first as it's a fundamental dependency for many packages
RUN pip install $(grep "^six" requirements.txt)

# install numpy first to avoid issues with bio python and bx-python (see also https://github.com/LUMC/vep2lovd/issues/1)
RUN pip install $(grep "^numpy" requirements.txt)

RUN pip install $(grep "^cython" requirements.txt)
RUN pip install -r requirements.txt -r test-requirements.txt

# install vcf tools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
RUN tar zxf vcftools*.tar.gz
RUN cd vcftools* && ./configure && make && make install && cd .. && rm -r vcftools*

# install BCFTools
RUN apt-get install -y bcftools
#RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2
#RUN tar -xjvf bcftools.tar.bz2
#RUN cd bcftools* && make && make prefix=/usr/local/bin install && cd .. && rm -r bcftools

# install tabix
RUN wget https://downloads.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2
RUN tar jxf tabix*.tar.bz2
RUN cd tabix* && make && cp tabix /usr/local/bin && cd .. && rm -r tabix*

ARG res=/files/resources

ENV BRCA_RESOURCES=$res

RUN mkdir -p $res /files/data && chmod -R o+rwx /files

RUN rm -r  /root/.cache

ARG FORCE_REBUILD=0
COPY . /opt/brca-exchange-kb

ENV LUIGI_CONFIG_PATH="/opt/luigi_pipeline_credentials.cfg"

ARG IS_GIT_DIRTY="False"
ARG GIT_COMMIT=""
LABEL GitCommit=$GIT_COMMIT IsGitDirty=$IS_GIT_DIRTY

CMD ["/opt/brca-exchange-kb/pipeline/docker/run_luigi.sh"]

