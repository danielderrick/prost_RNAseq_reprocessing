FROM ubuntu:17.04

ENV CUTADAPT_VERSION 1.10 
ENV TRIMGALORE_VERSION 0.4.3
ENV FASTQC_VERSION 0.11.5

RUN apt-get update
RUN apt-get install -y build-essential cython make openjdk-8-jdk python2.7-dev unzip wget

WORKDIR /opt

# Install cutadapt
RUN wget https://github.com/marcelm/cutadapt/archive/v1.10.tar.gz -O cutadaptv${CUTADAPT_VERSION}.tar.gz && \
    tar -zxvf cutadaptv${CUTADAPT_VERSION}.tar.gz && \
    cd cutadapt-${CUTADAPT_VERSION} && \
    python setup.py install && \
    python setup.py build_ext -i

# Install fastqc
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip -O fastqc-${FASTQC_VERSION}.zip && \
    unzip fastqc-${FASTQC_VERSION}.zip && \ 
    chmod +x FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/bin/

# Install TrimGalore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/${TRIMGALORE_VERSION}.tar.gz -O trim_galore-${TRIMGALORE_VERSION}.tar.gz && \
    tar -zxvf trim_galore-${TRIMGALORE_VERSION}.tar.gz && \
    cp TrimGalore-${TRIMGALORE_VERSION}/trim_galore /usr/bin/
