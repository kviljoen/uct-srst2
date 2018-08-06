# uct-cbio/uct-yamp Installation

To start using the uct-cbio/uct-yamp, follow the steps below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)

## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
sudo mv nextflow /usr/local/bin
```

### For Univeristy of Cape Town users working on HPC (hex):
```
#From your home directory on hex install nextflow
curl -fsSL get.nextflow.io | bash

#Add the following to ~/.bashrc:
JAVA_HOME=/opt/exp_soft/java/jdk1.8.0_31/
JAVA_CMD=/opt/exp_soft/java/jdk1.8.0_31/bin/java

#Do not run nextflow from the headnode, it requires substantial memory to run java. Please therefore first start an interactive job as follows: 
qsub -I -q UCTlong -l nodes=1:series600:ppn=1 -d `pwd`
```

**You need NextFlow version >= 0.24 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `uct-cbio/uct-yamp` is specified as the pipeline name when executing `nextflow run uct-cbio/uct-yamp`. If for some reason you need to use the development branch, this can be specified as `nextflow run uct-cbio/uct-yamp -r dev`

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/uct-cbio/uct-yamp/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/uct-yamp
```
## 3) Other requirements

YAMP requires a set of databases that are queried during its execution. Some of them should be automatically downloaded when installing the tools listed in the dependencies (or using specialised scripts, as those available with HUMAnN2), whilst other should be created by the user. Specifically, you will need:

- a FASTA file listing the adapter sequences to remove in the trimming step. This file should be available within the BBmap installation. If not, please download it from [here](https://github.com/BioInfoTools/BBMap/blob/master/resources/adapters.fa);
- two FASTA file describing synthetic contaminants. These files (`sequencing_artifacts.fa.gz` and `phix174_ill.ref.fa.gz`) should be available within the BBmap installation. If not, please download them from [here](https://sourceforge.net/projects/bbmap/);
- a FASTA file describing the contaminating genome(s). This file should be created by the users according to the contaminants present in their dataset. When analysing human metagenome, we suggest the users to always include the human genome. Please note that this file should be indexed beforehand. This can be done using BBMap, using the following command: `bbmap.sh -Xmx24G ref=my_contaminants_genomes.fa.gz `. 
	We suggest to download the FASTA file provided by Brian Bushnell for removing human contamination, using the instruction available [here](http://seqanswers.com/forums/showthread.php?t=42552);
- the BowTie2 database file for MetaPhlAn2. This file should be available within the MetaPhlAn2 installation. If not, please download it from [here](https://bitbucket.org/biobakery/metaphlan2/src/40d1bf693089836b5895623dd9ab1b21eb9a794c/db_v20/);
- the ChocoPhlAn and UniRef databases, that can be downloaded directly by HUMAnN2, as explained [here](https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-5-download-the-databases);
- [optional] a phylogenetic tree used by QIIME to compute a set of alpha-diversity measures (see [here](http://qiime.org/scripts/alpha_diversity.html) for details).

You can find an example of the folders layouts in this [wiki](https://github.com/alesssia/YAMP/wiki/Folders-layout-example) page.

You can also download all these files (please note that it might be necessary to edit this file list according to the analysis at hand) either from Zenodo ([https://zenodo.org/record/1068229#.Wh7a3rTQqL4](https://zenodo.org/record/1068229#.Wh7a3rTQqL4)), or using the following command:

```
wget https://zenodo.org/record/1068229/files/YAMP_resources_20171128.tar.gz
```

If you use this data file, please note that, before running YAMP, the FASTA file describing the human (contaminating) genome should be indexed with the following command:

```
bbmap.sh -Xmx24G ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
```

**Please also note that the size of this compressed data file is 16.7 GB.** 

## 4) Docker and/or Singularity setup
If you are not working on UCT hex, you may have to adapt the dockerfile in this repository for your own use, e.g. to add user-defined bind points. The `Dockerfile` can be built into a docker image on your own system (where docker has been installed) as follows:

First pull the git repository:
```
git clone https://github.com/uct-cbio/uct-yamp.git
```

Now build a local image by navigating to the folder where the `Dockerfile` is located, then run the following command (be careful to add the dot!):

```
docker build -t yampdocker .
```
If you are working on a cluster environment you will likely have to convert the docker image to a singularity image. This can be done as follows:

```
docker run -v /var/run/docker.sock:/var/run/docker.sock -v /home/katie/h3abionet16S/singularity-containers/:/output --privileged -t --rm singularityware/docker2singularity d02667d8d22e
```
Where `/home/katie/h3abionet16S/singularity-containers/` is the location where you want to save your singularity image and `d02667d8d22e` is the docker image ID obtained via `docker images` command

Next, test the singularity image:

```
singularity exec /scratch/DB/bio/singularity-containers/d02667d8d22e-2018-07-23-251e39cb1b13.img /bin/bash
```

Where `d02667d8d22e-2018-07-23-251e39cb1b13.img` is your singularity image. You are now in the singularity image environment and can test whether all software was successfully installed e.g. humann2 --help should print the relevant helpfile.

## 5) Optional manual setup
If you do not want to use the (recommended) containerised version of YAMP, you will need to install several additional tools for YAMP to work properly, and all of them should either be in the system path with execute and read permission, or made available within a multi-image scenario as the one we describe in the [multi-image scenario tutorial](https://github.com/alesssia/YAMP/wiki/multi-image-scenario).

The list of tools that should be available includes:
- fastQC v0.11.2+ ([http://www.bioinformatics.babraham.ac.uk/projects/fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc))
- BBmap v36.92+ ([https://sourceforge.net/projects/bbmap](http://www.bioinformatics.babraham.ac.uk/projects/fastqc))
- Samtools v1.3.1 ([http://samtools.sourceforge.net](http://samtools.sourceforge.net))
- MetaPhlAn2 v2.0+ ([https://bitbucket.org/biobakery/metaphlan2](https://bitbucket.org/biobakery/metaphlan2))
- QIIME v1.9.1+ ([http://qiime.org](http://qiime.org))
- HUMAnN2 v0.9.9+ ([https://bitbucket.org/biobakery/humann2](https://bitbucket.org/biobakery/humann2))

Following the links, you will find detailed instructions on how to install them, as explained by their developers. 
Notably, MetaPhlAn2, QIIME, and HUMAnN2 are also available in [bioconda](https://anaconda.org/bioconda/). 

---

[![UCT Computational Biology](/assets/cbio_logo.png)](http://www.cbio.uct.ac.za/)

---
