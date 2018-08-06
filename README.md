# ![kviljoen/YAMP](/assets/cbio_logo.png)

# UCT-YAMP: A reference-based WGS metagenomics pipeline, adapted from Yet Another Metagenomic Pipeline (YAMP), implemented in Nextflow

A reference-based WGS metagenomics pipeline using the Nextflow workflow manager. This pipeline accepts raw reads in .fastq format, performs quality filtering, adapter removal and decontamination, followed by taxonomic profiling with MetaPhlAn2, and functional profiling with HUMAnN2. This pipeline was adapted from https://github.com/alesssia/YAMP for implementation on the University of Cape Town (UCT) high-performance compute cluster.

## Basic usage:

    The typical command for running the pipeline is as follows:
    nextflow run uct-cbio/uct-yamp --reads '*_R{1,2}.fastq.gz' -profile uct_hex
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. uct_hex OR standard
    Options:
      --singleEnd                   Specifies that the input is single end reads
    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
    Trimming options
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run 					  sent to you when the workflow exits
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.  
     Help:
      --help                        Will print out summary above when executing nextflow run uct-cbio/16S-rDNA-dada2-pipeline                                     --help

## Prerequisites

Nextflow (0.26.x or higher), all other software/tools required are contained in the (platform-independent) dockerfile, which is converted to a singularity image for use on a cluster environment.

## Documentation
The uct-cbio/uct-yamp pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)

## Built With

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/what-docker)
* [Singularity](https://singularity.lbl.gov/)


## Credits

The YAMP pipeline was built by Dr Alessia Visconti (https://github.com/alesssia/YAMP). Please remember to cite Dr Visconti  
> Visconti A,. Martin T.C., and Falchi M., *"YAMP: a containerised workflow enabling reproducibility in metagenomics research"*, GigaScience (2018), [https://doi.org/10.1093/gigascience/giy072](https://doi.org/10.1093/gigascience/giy072)
when using this pipeline. Further development to the Nextflow workflow and containerisation in Docker and Singularity for implementation specifically on UCT's HPC was done by Dr Katie Lennard, with Nextflow template inspiration and code snippets from Phil Ewels http://nf-co.re/

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


