FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    git \
    liblzo2-dev \
    pkg-config \
    python-pip \
    python-gdbm \
    python-lzo \
    vim \
    wget \
    zlib1g-dev

WORKDIR /opt

COPY requirements_docker.txt .

# install numpy first to avoid issues with bio python and bx-python (see also https://github.com/LUMC/vep2lovd/issues/1)
RUN pip install $(grep numpy requirements_docker.txt)

RUN pip install -r requirements_docker.txt

# install vcf tools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz
RUN tar zxf vcftools*.tar.gz
RUN cd vcftools* && ./configure && make && make install && cd .. && rm -r vcftools*

# install tabix
RUN wget https://downloads.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2
RUN tar jxf tabix*.tar.bz2
RUN cd tabix* && make && cp tabix /usr/local/bin && cd .. && rm -r tabix*

ARG res=/files/resources

ENV BRCA_RESOURCES=$res

RUN mkdir -p $res /files/data && chmod -R o+rwx /files

RUN mkdir -p /.synapseCache && chmod o+rwx /.synapseCache  

RUN git clone https://github.com/counsyl/hgvs.git
# taking pyhgvs 0.9.4
RUN cd hgvs && git checkout aebe5bd9683f4b5937fd653ce4b13fcd4f3ebb10 && python setup.py install


#leiden brca
RUN git clone https://github.com/BRCAChallenge/leiden.git && cd leiden && python setup.py install

RUN rm -r /opt/leiden /opt/hgvs /root/.cache

ARG FORCE_REBUILD=0
ARG BRCA_GIT_REPO=https://github.com/BRCAChallenge/brca-exchange.git
ARG BRCA_EXCHANGE_COMMIT=master
RUN git clone ${BRCA_GIT_REPO} && cd brca-exchange && git checkout ${BRCA_EXCHANGE_COMMIT}

COPY run_luigi.sh .

ENV LUIGI_CONFIG_PATH="/opt/luigi_pipeline_credentials.cfg"

CMD ["/opt/run_luigi.sh"]

