# ![kviljoen/YAMP](/assets/cbio_logo.png)
# uct-yamp usage

## General Nextflow info
Nextflow handles job submissions on different environments (PBS in our case), and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` or a similar tool. 

On UCT's hex cluster you would start a screen session from the headnode and then start and interactive job. Once you are on a worker node you can enter the typical command for running the pipeline specified below.

It is recommended that you limit the Nextflow Java virtual machine's memory usage by adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
For University of Cape Town users who will be running Nextflow on UCT's HPC (hex), you need to include the following two lines in your `~/.bashrc`:

```bash
JAVA_CMD=/opt/exp_soft/java/jdk1.8.0_31/bin/java
export PATH=$PATH:/opt/exp_soft/cbio/nextflow
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run uct-cbio/uct-yamp --reads '*_R{1,2}.fastq.gz' -profile uct_hex
```

This will launch the pipeline with the `uct_hex` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
outdir         # Final results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull uct-cbio/uct-yamp
```
Note that `nextflow pull` by default stores pipelines in $HOME/.nextflow/assets To clone this pipeline to a specific directory, use `nextflow clone uct-cbio/uct-yamp target-dir ` where target-dir is the target directory.

## Main (required) arguments

### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different compute environment - follow the links below to see instructions for running on that system. Available profiles are:

* `standard`
    * The default profile, used if `-profile` is not specified at all. Runs locally and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles.
* `uct_hex`
    * Designed to run on UCT's high-performance cluster (hex).
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).
    
* You can build your own profile suited to your cluster. See [`conf/uct_hex.config`](../conf/uct_hex.config) as an example

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_R{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

## Paths to external resources
### Adapter sequences and synthetic contaminants to be removed in the trimming step
Adapters are removed using the BBTools package, currently includes Illumina Truseq and Nextera adapters sequences in our singularity installation of /opt/conda/opt/bbmap-37.10/resources/adapters.fa. You can specify whether or not BBDuk looks for the reverse-complement of the reference sequences as well as the forward sequence with the flag “rcomp=t” or “rcomp=f”; by default it looks for both. You can also specify custom adapters if necessary (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
Also available in the /opt/conda/opt/bbmap-37.10/resources/ folder are sequencing_artifacts.fa.gz and phix174_ill.ref.fa.gz 

### Reference pan-genome for decontamination
A FASTA file describing the contaminating genome(s). This file should be created according to the contaminants present in your dataset. When analysing the human metagenome, include the human genome and specify the file path as `refForeingGenome` in the config file. Please note that this file should be indexed beforehand. This can be done using BBMap, using the following command: bbmap.sh -Xmx24G ref=my_contaminants_genomes.fa.gz. For the default UCT option (see conf/uct_hex.config) the human contaminants file was downloaded from https://zenodo.org/record/1208052/files/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz

### BowTie2 database for MetaPhlAn2


## Job resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for an example.

## General command line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the summary HTML / e-mail.

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, if you don't want FastQC errors to be ignored, you can specify a config file using `-c` that contains the following:

```groovy
process.$fastqc.errorStrategy = 'terminate'
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.



