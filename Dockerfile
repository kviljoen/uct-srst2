FROM debian:jessie

ENV CONDA_INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
#Exports conda path
ENV PATH $PATH:/opt/conda/bin/

#Installs packages that will be used either for building other software 
#or directly by YAMP
#RUN apk --update add --no-cache bash procps wget curl gzip perl mesa-gl 

#Installs miniconda 

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    wget \
    ca-certificates \
    bash \
    procps \
    wget \
    curl \
    gzip \
    perl \
    mesa-gl && \
    wget --quiet https://repo.continuum.io/miniconda/${CONDA_INSTALLER} && \
    /bin/bash /${CONDA_INSTALLER} -b -p /opt/conda && \
    rm ${CONDA_INSTALLER} && \
    /opt/conda/bin/conda install --yes conda && \
    conda install conda-build && \
    conda remove tk --yes && \
    conda clean --yes --tarballs --packages --source-cache && \
    apt-get purge -y --auto-remove wget ca-certificates  && \
    apt-get clean


#Update conda and uses it to install software used by YAMP
#that is required to use YAMP on AWS Batch
RUN conda install -c bioconda -y bbmap=37.10 fastqc=0.11.5 metaphlan2=2.6.0 qiime=1.9.1 humann2=0.9.9 &&\
    conda install -c conda-forge -y awscli &&\







