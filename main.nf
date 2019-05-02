#!/usr/bin/env nextflow
/*
========================================================================================
              		   UCT-SRST2   P I P E L I N E
========================================================================================
 			MLST NEXTFLOW PIPELINE USING SRST2
 
----------------------------------------------------------------------------------------
*/

/**
	Prints help when asked for
*/

def helpMessage() {
    log.info"""
    ===================================
     uct-srst2  ~  version ${params.version}
    ===================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run uct-cbio/uct-yamp --reads '*_R{1,2}.fastq.gz' -profile uct_hex
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. uct_hex OR standard
      --mlst_db			    Specify fasta-formatted mlst database for srst2
      --mlst_definitions	    Definitions for MLST scheme (required if mlst_db
                        	    supplied and you want to calculate STs)
    Other srst2 options:
      --mlst_delimiter		    Character(s) separating gene name from allele number
                        	    in MLST database (default "-", as in arcc-1)
      --mlst_max_mismatch	    Maximum number of mismatches per read for MLST allele calling (default 10)
      --AMR_db			    Fasta-fromatted gene databases for resistance gene analysis (optional)
      
    General options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      
     Help:
      --help                        Will print out summary above when executing nextflow run uct-cbio/uct-yamp --help                                    
    """.stripIndent()
}
	
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
//params.project = false
params.email = false
params.plaintext_email = false

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//Validate inputs	

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Validate user-inputs
if( params.mlst_db ) {
        mlst_db = file(params.mlst_db)
        if( !mlst_db.exists() ) exit 1, "MLST DB file could not be found: ${params.mlst_db}"
}

if( params.mlst_definitions ) {
        mlst_definitions= file(params.mlst_definitions)
        if( !mlst_definitions.exists() ) exit 1, "MLST definitions file could not be found: ${params.mlst_definitions}"
}

// Returns a tuple of read pairs in the form
// [sample_id, forward.fq, reverse.fq] where
// the dataset_id is the shared prefix from
// the two paired FASTQ files.
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }
    .into { ReadPairsToSrst2 }

// Header log info
log.info "==================================="
log.info " uct-srst2  ~  version ${params.version}"
log.info "==================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['OS']		= System.getProperty("os.name")
summary['OS.arch']	= System.getProperty("os.arch") 
summary['OS.version']	= System.getProperty("os.version")
summary['javaversion'] = System.getProperty("java.version") //Java Runtime Environment version
summary['javaVMname'] = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
summary['javaVMVersion'] = System.getProperty("java.vm.version") //Java Virtual Machine implementation version
//Gets starting time		
sysdate = new java.util.Date() 
summary['User']		= System.getProperty("user.name") //User's account name
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

		
/*
 *
 * Step 1: srst2 (run per sample)
 *
 */

process srst2 {
    tag { "srst2.${pairId}" }
    publishDir "${params.outdir}/srst2", mode: "copy", overwrite: false

    input:
        set pairId, file(reads) from ReadPairsToSrst2

    output:
        file "${pairId}_srst2*results.txt" into srst2_results
	
    script:
    
    """
    srst2 --forward ${reads[0]} --reverse ${reads[1]} --output ${pairId}_srst2 --mlst_db $mlst_db \
    --mlst_definitions $mlst_definitions --mlst_delimiter $params.mlst_delimiter --gene_db $AMR_db
    """
}





