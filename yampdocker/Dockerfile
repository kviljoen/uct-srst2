FROM continuumio/miniconda
MAINTAINER Katie Lennard
LABEL authors="katieviljoen@gmail.com" \
    description="Docker image containing all requirements for the uct-cbio/YAMP pipeline"

RUN apk --update add --no-cache bash procps wget curl gzip perl mesa-gl 

COPY environment.yml /
RUN conda update -n base conda && \
    conda env create -f /environment.yml && \
    conda clean -a

ENV PATH /opt/conda/envs/YAMP/bin:$PATH
################## UCT Hex specific ###########################
RUN mkdir -p /researchdata/fhgfs
RUN mkdir -p /scratch/DB/bio/YAMP
