# ![kviljoen/YAMP](/assets/cbio_logo.png)

# UCT-YAMP: A reference-based WGS metagenomics pipeline, adapted from Yet Another Metagenomic Pipeline (YAMP), implemented in Nextflow

A reference-based WGS metagenomics pipeline using the Nextflow workflow manager. This pipeline accepts raw reads in .fastq format, performs quality filtering, adapter removal and decontamination, followed by taxonomic profiling with MetaPhlAn2, and functional profiling with HUMAnN2. This pipeline was adapted from https://github.com/alesssia/YAMP for implementation on the University of Cape Town (UCT) high-performance compute cluster.

## Basic usage:

    The typical command for running the pipeline is as follows:
    nextflow run uct-cbio/uct-yamp --reads '*_R{1,2}.fastq.gz' -profile uct_hex

    Mandatory arguments:
      --reads			Path to input data (must be surrounded with quotes)
      -profile			Hardware config to use. uct_hex OR standard
      
    BBduk trimming options:
      --qin			Input quality offset: 33 (ASCII+33) or 64 (ASCII+64, default=33
      --kcontaminants		Kmer length used for finding contaminants, default=23	
      --phred			Regions with average quality BELOW this will be trimmed, default=10 
      --minlength		Reads shorter than this after trimming will be discarded, default=60
      --mink			Shorter kmers at read tips to look for, default=11 
      --hdist			Maximum Hamming distance for ref kmers, default=1            
    BBwrap parameters for decontamination:	
      --mind			Approximate minimum alignment identity to look for, default=0.95
      --maxindel		Longest indel to look for, default=3
      --bwr			Restrict alignment band to this, default=0.16
	
    MetaPhlAn2 parameters: 
      --bt2options		Presets options for BowTie2, default="very-sensitive"
      
    Other options:
      --keepQCtmpfile		Whether the temporary files resulting from QC steps should be kept, default=false
      --keepCCtmpfile		Whether the temporary files resulting from MetaPhlAn2 and HUMAnN2 should be kept, default=false 
      --outdir			The output directory where the results will be saved
      --email                   Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      
     Help:
      --help			Will print out summary above when executing nextflow run uct-cbio/uct-yamp --help 

## Prerequisites

Nextflow (0.26.x or higher), all other software/tools required are contained in the (platform-independent) dockerfile, which is converted to a singularity image for use on a cluster environment.

There are however reference databases that need to be downloaded (only if you are NOT working on UCT hex), these are:
1. Reference pan-genome for contamination (used by process decontaminate) can be downloaded from here: https://zenodo.org/record/1208052/files/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
2. BowTie2 database for MetaPhlAn2: The mpa_v20_m200.pkl can be downloaded from here https://bitbucket.org/biobakery/metaphlan2/downloads/ (mpa_v20_m200.tar file); the bowtie2 index files for metaphlan2 (which should be specified under "bowtie2db=" in the uct_hex.config file can be downloaded from here https://www.dropbox.com/sh/w4j4yr1b0o7xu9v/AAAx1yiV6enIGR7SuC8B34cKa?dl=0 Note that these index files can also be built from the fna file in the .tar in step 1. as outlined here https://groups.google.com/forum/#!topic/metaphlan-users/5ltD8X8Xitc

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


