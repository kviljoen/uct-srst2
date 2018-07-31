FROM continuumio/miniconda
MAINTAINER Katie Lennard
LABEL authors="katieviljoen@gmail.com" \
    description="Docker image containing all requirements for the uct-cbio/YAMP pipeline"

#install basic unix stuff
RUN apt-get update && apt-get install -y \
bash procps wget curl gzip perl \
&& rm -rf /var/lib/apt/lists/*

#install specialized pipeline packages
COPY environment.yml /
RUN conda update -n base conda && \
    conda env create -f /environment.yml && \
    conda clean -a

ENV PATH /opt/conda/envs/YAMP/bin:$PATH

################## UCT Hex specific ###########################
RUN mkdir -p /researchdata/fhgfs
RUN mkdir -p /scratch/DB/bio/YAMP
