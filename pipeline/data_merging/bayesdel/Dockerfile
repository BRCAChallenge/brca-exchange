FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    curl \
    parallel
    #&& rm -rf /var/lib/apt/lists/*

RUN mkdir /opt/victor

WORKDIR /opt/victor

RUN curl -L -o VICTOR_linux.tgz https://www.dropbox.com/s/rkfu47k8o4mh5wg/VICTOR_linux.tgz && tar zxf VICTOR_linux.tgz && rm -rf VICTOR_linux.tgz

# TODO split up and do it manually

# installing htslib (includes tabix and bgzip)

WORKDIR /opt

ARG HTSLIB_VERSION=1.10.2
# TODO squash it together
RUN curl -L -o htslib-$HTSLIB_VERSION.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2

RUN tar -vxjf htslib-$HTSLIB_VERSION.tar.bz2

RUN apt-get install -y \
    make \
    g++ \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-nss-dev

RUN cd htslib-$HTSLIB_VERSION && make

RUN cd htslib-$HTSLIB_VERSION && make install

# tabix bgzip parallel
#RUN curl -sL http://www.bjfenglab.org/download/get_pkg_inst.sh | bash
ENV PATH=/opt/victor/VICTOR:${PATH}
# setup depencencies
# curl -sL http://www.bjfenglab.org/download/get_pkg_inst.sh | bash
# tweak env variables

COPY run_annotation.sh .
RUN chmod +rx run_annotation.sh

RUN mkdir /opt/victor/wdir && chmod o+rwx /opt/victor/wdir
WORKDIR /opt/victor/wdir

USER nobody

ENTRYPOINT ["/opt/run_annotation.sh"]

# --> data depencency
