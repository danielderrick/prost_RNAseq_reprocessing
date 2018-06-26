FROM ubuntu:16.04

ENV KALLISTO_VERSION 0.43.0

RUN apt-get update
RUN apt-get install -y build-essential cmake gcc hdf5-tools libhdf5-dev make wget zlib1g-dev

WORKDIR /opt

RUN wget https://github.com/pachterlab/kallisto/archive/v${KALLISTO_VERSION}.tar.gz -O /opt/kallisto-v${KALLISTO_VERSION}.tar.gz && \ 
    tar -zxvf kallisto-v${KALLISTO_VERSION}.tar.gz && \
    cd kallisto-${KALLISTO_VERSION} && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install
