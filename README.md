# ![kviljoen/uct-srst2](/assets/cbio_logo.png)

# UCT-srst2: MLST (and antimicrobial resistance/ virulence factor) analysis using srst2, implemented in Nextflow

This pipeline accepts raw (or trimmed) reads in .fastq format and performs MLST (and optional AMR or virulence factor) analyses using srst2 (https://github.com/katholt/srst2). Implemented for the Ilifu (slurm) cluster, which can be modified by adding a custom config file. The srst2 singularity image used was obtained here https://quay.io/repository/biocontainers/srst2?tab=tags

## Basic usage:

    The typical command for running the pipeline is as follows:
    nextflow run kviljoen/uct-srst2 --reads '*_R{1,2}.fastq.gz' -profile ilifu --mlst_db Pseudomonas_aeruginosa.fasta --mlst_definitions paeruginosa.txt --mlst_delimiter '_' --AMR_db ARGannot_r3.fasta --outdir Ps_aerug_srst2_MLST

    Mandatory arguments:
      --reads			Path to input data (must be surrounded with quotes!)
      -profile			Hardware config to use. ilifu OR standard
      --outdir			Specify folder where results should be written to (will be created if non-existent)
    
    srst2 options:
      --mlst_delimiter		Default="-" you may have to change to e.g. "_" depending on the format of your mlst reference file
      --mlst_db			Fasta file of MLST alleles (can be downloaded with getmlst.py from srst2 e.g. getmlst.py --species "Escherichia coli"   
      --gene_db			Antimicrobial resistence (or other e.g virulence) gene DB (Fasta) (can be downloaded from https://github.com/katholt/srst2/tree/master/data or created)
      --min_gene_cov    Minimum %coverage cutoff for gene reporting (default 90)
      --max_gene_divergence Maximum %divergence cutoff for gene reporting (default 10)
      

     Help:
      --help			Will print out summary above when executing nextflow run kviljoen/uct-srst2 --help 

## Prerequisites

Nextflow (0.27.6 or higher), all other software/tools required are contained in the (platform-independent) singularity image for srst2, which can be obtained here https://quay.io/repository/biocontainers/srst2?tab=tags


## Documentation
For detailed information please see https://github.com/katholt/srst2

## Built With

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/what-docker)
* [Singularity](https://singularity.lbl.gov/)


## Credits
Please remember to cite the authors of srst2, see here http://genomemedicine.com/content/6/11/90

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


