#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
               M E T A G E N O M I C S   P I P E L I N E
========================================================================================
 METAGENOMICS NEXTFLOW PIPELINE ADAPTED FROM YAMP FOR UCT CBIO
 
----------------------------------------------------------------------------------------
*/

/**
	Prints help when asked for
*/

def helpMessage() {
    log.info"""
    ===================================
     uct-cbio/uct-yamp  ~  version ${params.version}
    ===================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run uct-cbio/uct-yamp --reads '*_R{1,2}.fastq.gz' -profile uct_hex
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. uct_hex OR standard
      
    BBduk trimming options:
      --qin			    Input quality offset: 33 (ASCII+33) or 64 (ASCII+64, default=33
      --kcontaminants		    Kmer length used for finding contaminants, default=23	
      --phred			    Regions with average quality BELOW this will be trimmed, default=10 
      --minlength		    Reads shorter than this after trimming will be discarded, default=60
      --mink			    Shorter kmers at read tips to look for, default=11 
      --hdist			    Maximum Hamming distance for ref kmers, default=1            

    BBwrap parameters for decontamination:	
      --mind			   Approximate minimum alignment identity to look for, default=0.95
      --maxindel		   Longest indel to look for, default=3
      --bwr			   Restrict alignment band to this, default=0.16
	
    MetaPhlAn2 parameters: 
      --bt2options 		   Presets options for BowTie2, default="very-sensitive"
      
    Other options:
      --keepQCtmpfile		    Whether the temporary files resulting from QC steps should be kept, default=false
      --keepCCtmpfile		    Whether the temporary files resulting from MetaPhlAn2 and HUMAnN2 should be kept, default=false 
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

//if (params.librarylayout != "paired" && params.librarylayout != "single") { 
//	exit 1, "Library layout not available. Choose any of <single, paired>" 
//}   

if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { ReadPairsToQual; ReadPairs }

// Header log info
log.info "==================================="
log.info " uct-cbio/uct-yamp  ~  version ${params.version}"
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
 * Step 1: FastQC (run per sample)
 *
 */

process runFastQC {
    tag { "rFQC.${pairId}" }
    publishDir "${params.outdir}/FilterAndTrim", mode: "copy", overwrite: false

    input:
        set pairId, file(in_fastq) from ReadPairsToQual

    output:
        file("${pairId}_fastqc/*.zip") into fastqc_files

    """
    mkdir ${pairId}_fastqc
    fastqc --outdir ${pairId}_fastqc \
    ${in_fastq.get(0)} \
    ${in_fastq.get(1)}
    """
}

process runMultiQC{
    tag { "rMQC" }
    publishDir "${params.outdir}/FilterAndTrim", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

/*
 *
 * Step 2: De-duplication (run per sample)
 *
 */

process dedup {
	tag{ "dedup" }

	input:
	set val(pairId), file(reads) from ReadPairs

	output:
	set val(pairId), file("${pairId}_dedupe_R1.fq"), file("${pairId}_dedupe_R2.fq") into totrim, topublishdedupe

	script:
	"""
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	clumpify.sh -Xmx\"\$maxmem\" in1="${reads[0]}" in2="${reads[1]}" out1=${pairId}_dedupe_R1.fq out2=${pairId}_dedupe_R2.fq \
	qin=$params.qin dedupe subs=0 threads=${task.cpus}
	
	"""
}


/*
 *
 * Step 3: BBDUK: trim + filter (run per sample)
 *
 */

process bbduk {
	tag{ "bbduk" }
	
	input:
	set val(pairId), file("${pairId}_dedupe_R1.fq"), file("${pairId}_dedupe_R2.fq") from totrim
	file(adapters) from Channel.from( file(params.adapters) )
	file(artifacts) from Channel.from( file(params.artifacts) )
	file(phix174ill) from Channel.from( file(params.phix174ill) )

	output:
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq"), file("${pairId}_trimmed_singletons.fq") into todecontaminate, topublishtrim
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq") into filteredReadsforQC

	script:
	"""	
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	#Quality and adapter trim:
	bbduk.sh -Xmx\"\$maxmem\" in=${pairId}_dedupe_R1.fq in2=${pairId}_dedupe_R2.fq out=${pairId}_trimmed_R1_tmp.fq \
	out2=${pairId}_trimmed_R2_tmp.fq outs=${pairId}_trimmed_singletons_tmp.fq ktrim=r \
	k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred \
	minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe 
	
	#Synthetic contaminants trim:
	bbduk.sh -Xmx\"\$maxmem\" in=${pairId}_trimmed_R1_tmp.fq in2=${pairId}_trimmed_R2_tmp.fq \
	out=${pairId}_trimmed_R1.fq out2=${pairId}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts \
	qin=$params.qin threads=${task.cpus} 

	#Synthetic contaminants trim for singleton reads:
	bbduk.sh -Xmx\"\$maxmem\" in=${pairId}_trimmed_singletons_tmp.fq out=${pairId}_trimmed_singletons.fq \
	k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=${task.cpus}

	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${pairId}_trimmed*_tmp.fq 

	"""
}


/*
 *
 * Step 4: FastQC post-filter and -trim (run per sample)
 *
 */

process runFastQC_postfilterandtrim {
    tag { "rFQC_post_FT.${pairId}" }
    publishDir "${params.outdir}/FastQC_post_filter_trim", mode: "copy", overwrite: false

    input:
    	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq") from filteredReadsforQC

    output:
        file("${pairId}_fastqc_postfiltertrim/*.zip") into fastqc_files_2

    """
    mkdir ${pairId}_fastqc_postfiltertrim
    fastqc --outdir ${pairId}_fastqc_postfiltertrim \
    ${pairId}_trimmed_R1.fq \
    ${pairId}_trimmed_R2.fq
    """
}

process runMultiQC_postfilterandtrim {
    tag { "rMQC_post_FT" }
    publishDir "${params.outdir}/FastQC_post_filter_trim", mode: 'copy', overwrite: false

    input:
        file('*') from fastqc_files_2.collect()

    output:
        file('multiqc_report.html')

    """
    multiqc .
    """
}

/*
 *
 * Step 5: Decontamination (run per sample)
 *
 */

process decontaminate {
	tag{ "decon" }
	publishDir  "${params.outdir}/decontaminate" , mode: 'move', pattern: "*_clean.fq.gz"
	
	input:
	set val(pairId), file("${pairId}_trimmed_R1.fq"), file("${pairId}_trimmed_R2.fq"), file("${pairId}_trimmed_singletons.fq") from todecontaminate
	file(refForeingGenome) from Channel.from( file(params.refForeingGenome, type: 'dir') )
	
	output:
	file "*_clean.fq.gz"
	set val(pairId), file("${pairId}_clean.fq") into cleanreadstometaphlan2, cleanreadstohumann2 
	set val(pairId), file("${pairId}_cont.fq") into topublishdecontaminate
	
	script:
	"""
	
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	#Decontaminate from foreign genomes
	bbwrap.sh  -Xmx\"\$maxmem\" mapper=bbmap append=t in1=${pairId}_trimmed_R1.fq,${pairId}_trimmed_singletons.fq in2=${pairId}_trimmed_R2.fq,null \
	outu=${pairId}_clean.fq outm=${pairId}_cont.fq minid=$params.mind \
	maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred \
	path=$refForeingGenome qin=$params.qin threads=${task.cpus} untrim quickmatch fast
	
	gzip -c ${pairId}_clean.fq > ${pairId}_clean.fq.gz

	"""
}


/*
 *
 * Step 6:  metaphlan2 (run per sample)
 *
 */

process metaphlan2 {
	tag{ "metaphlan2" }
	
	publishDir  "${params.outdir}/metaphlan2", mode: 'copy', pattern: "*.{biom,tsv}", overwrite: false
	
	input:
	set val(pairId), file(infile) from cleanreadstometaphlan2
	file(mpa_pkl) from Channel.from( file(params.mpa_pkl) )
	file(bowtie2db) from Channel.fromPath( params.bowtie2db, type: 'dir' )

    output:
	file "${pairId}_metaphlan_profile.tsv" into metaphlantohumann2
	file "${pairId}_bt2out.txt" into topublishprofiletaxa


	script:
	"""
	#If a file with the same name is already present, Metaphlan2 will crash
	rm -rf ${pairId}_bt2out.txt
	
	#Estimate taxon abundances
	metaphlan2.py --input_type fastq --tmp_dir=. --biom ${pairId}.biom --bowtie2out=${pairId}_bt2out.txt \
	--mpa_pkl $mpa_pkl  --bowtie2db $bowtie2db/$params.bowtie2dbfiles --bt2_ps $params.bt2options --nproc ${task.cpus} \
	$infile ${pairId}_metaphlan_profile.tsv
	
	#Sets the prefix in the biom file
	sed -i 's/Metaphlan2_Analysis/${pairId}/g' ${pairId}.biom
	sed -i 's/Metaphlan2_Analysis/${pairId}/g' ${pairId}_metaphlan_profile.tsv
	
	#KL add: make one file combining all samples (needs testing, perphaps own process)
	merge_metaphlan_tables.py *_metaphlan_profile.tsv > metaphlan_merged_abundance_table.tsv

	"""
}


process humann2 {

	publishDir  "${params.outdir}/humann2", mode: 'copy', pattern: "*.{tsv,log}", overwrite: false
	
	input:
	set val(pairId), file(cleanreads) from cleanreadstohumann2
	file(humann2_profile) from metaphlantohumann2
	file(chocophlan) from Channel.fromPath( params.chocophlan, type: 'dir' )
	file(uniref) from Channel.fromPath( params.uniref, type: 'dir' )
	
    output:
	file "${pairId}_genefamilies.tsv"
	file "${pairId}_pathcoverage.tsv"
	file "${pairId}_pathabundance.tsv"
	
	//Those may or may not be kept, according to the value of the keepCCtmpfile parameter
	set ("${pairId}_bowtie2_aligned.sam", "${pairId}_bowtie2_aligned.tsv", "${pairId}_diamond_aligned.tsv", 
	     "${pairId}_bowtie2_unaligned.fa", "${pairId}_diamond_unaligned.fa") into topublishhumann2	

	script:
	"""
	#Functional annotation
	humann2 --input $cleanreads --output . --output-basename ${pairId} \
	--taxonomic-profile $humann2_profile --nucleotide-database $chocophlan --protein-database $uniref \
	--pathways metacyc --threads ${task.cpus} --memory-use minimum

	
	#Performs functional annotation, redirect is done here because HUMAnN2 freaks out


	#Some of temporary files (if they exist) may be moved in the working directory, 
	#according to the keepCCtmpfile parameter. Others (such as the bowties2 indexes), 
	#are always removed. Those that should be moved, but have not been created by 
	#HUMAnN2, are now created by the script (they are needed as output for the channel)
	files=(${params.prefix}_bowtie2_aligned.sam ${params.prefix}_bowtie2_aligned.tsv ${params.prefix}_diamond_aligned.tsv \
	${params.prefix}_bowtie2_unaligned.fa ${params.prefix}_diamond_unaligned.fa)
	
	for i in {1..5}
	do
		if [ -f ${params.prefix}_humann2_temp/\${files[((\$i-1))]} ]
		then
			mv ${params.prefix}_humann2_temp/\${files[((\$i-1))]} .
		else
			touch \${files[((\$i-1))]}
		fi
	done
	rm -rf ${params.prefix}_humann2_temp/

 	"""
}


/*
 *
 * Step 7:  Save trimmed and decontaminated reads if requested
 *
 */
	
	
process saveQCtmpfile {

	publishDir  "${params.outdir}/QCtmpfiles", mode: 'copy'
		
	input:
	file (tmpfile) from topublishdedupe.mix(topublishtrim, topublishdecontaminate).flatMap()

	output:
	file "*.fq.gz"

	when:
	params.keepQCtmpfile
		
	script:
	"""
	gzip --force -c $tmpfile > ${tmpfile}.gz
	"""
}

/*
 *
 * Step 8:  Save tmp files from metaphlan2 and humann2 if requested
 *
 */	
	
process saveCCtmpfile {

	publishDir  "${params.outdir}/CCtmpfiles", mode: 'copy'
		
	input:
	file (tmpfile) from topublishprofiletaxa.mix(topublishhumann2).flatMap()

	output:
	file "$tmpfile"

	when:
	params.keepCCtmpfile
		
	script:
	"""
	echo $tmpfile
	"""
}
