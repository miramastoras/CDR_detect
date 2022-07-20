FROM ubuntu:18.04
MAINTAINER Mira Mastoras mmastora@ucsc.edu

ARG git_commit

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev libhdf5-dev liblzma-dev && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev jq pigz \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 --no-cache-dir install --upgrade pip \
  && rm -rf /var/lib/apt/lists/*

RUN python -m pip install pandas pysam biopython numpy argparse

# get samtools
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
    rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
    cd samtools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses --disable-lzma && \
    make && \
    ln -s /opt/samtools/samtools-1.9/samtools /usr/local/bin/samtools

# set environment variables
ENV bedtools_version 2.29.2

# Install dependencies
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    vim \
    wget \
    curl \
    zlib1g-dev \
    python-pip \
    libbz2-dev \
    liblzma-dev

# install bedtools
WORKDIR /usr/local/bin
RUN curl -SL https://github.com/arq5x/bedtools2/archive/v${bedtools_version}.tar.gz \
    > v${bedtools_version}.tar.gz
RUN tar -xzvf v${bedtools_version}.tar.gz
WORKDIR /usr/local/bin/bedtools2-${bedtools_version}
RUN make
RUN ln -s /usr/local/bin/bedtools2-${bedtools_version}/bin/bedtools /usr/local/bin/bedtools

# install CDR detect
WORKDIR /opt
RUN git clone https://github.com/miramastoras/CDR_detect && \
    cd CDR_detect/scripts && \
    chmod +x CDR_detect.py &&
ENV PATH="/opt/CDR_detect/scripts/:${PATH}"

WORKDIR /data
