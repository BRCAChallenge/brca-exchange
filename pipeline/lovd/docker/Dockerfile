FROM ubuntu:18.04

RUN chmod 1777 /tmp

RUN apt-get update && apt-get install -y git \
    python-lzo \
    python-pip \
    liblzo2-dev \
    libbz2-dev \
    lzma-dev \
    liblzma-dev \
    zlib1g-dev

WORKDIR /opt

RUN pip install psycopg2-binary==2.8.4 configparser==3.5.0 backports.functools_lru_cache==1.6.1 bioutils==0.3.3 git+https://github.com/biocommons/biocommons.seqrepo.git@129d285727228920dd2d7746a47c48eb7b88b191
RUN git clone https://github.com/counsyl/hgvs.git
# taking pyhgvs 0.9.4
RUN cd hgvs && git checkout aebe5bd9683f4b5937fd653ce4b13fcd4f3ebb10 && python setup.py install

#leiden brca
RUN git clone https://github.com/BRCAChallenge/leiden.git && cd leiden && git checkout d5352801da0858840d87280f36dbce14159a6dd4 && python setup.py install

RUN rm -r /opt/leiden /root/.cache

RUN mkdir /data && chmod -R o+rwx /data

ENTRYPOINT ["extract_data.py"]

