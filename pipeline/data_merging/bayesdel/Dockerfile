FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    curl \
    g++ \
    libbz2-dev \
    libcurl4-nss-dev \
    liblzma-dev \
    make \
    parallel \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir /opt/victor

WORKDIR /opt/victor

RUN curl -L -o VICTOR_linux.tgz https://www.dropbox.com/s/rkfu47k8o4mh5wg/VICTOR_linux.tgz && tar zxf VICTOR_linux.tgz && rm -rf VICTOR_linux.tgz

WORKDIR /opt

ARG HTSLIB_VERSION=1.10.2

# installing htslib (includes tabix and bgzip)
RUN curl -L -o htslib-$HTSLIB_VERSION.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 \
  && tar -vxjf htslib-$HTSLIB_VERSION.tar.bz2 \
  && cd htslib-$HTSLIB_VERSION \
  && make \
  && make install

ENV PATH=/opt/victor/VICTOR:${PATH}

COPY run_annotation.sh .

RUN chmod +rx run_annotation.sh

RUN mkdir /opt/victor/wdir && chmod o+rwx /opt/victor/wdir
WORKDIR /opt/victor/wdir

USER nobody

ENTRYPOINT ["/opt/run_annotation.sh"]
