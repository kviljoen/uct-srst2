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
    Options:
      --singleEnd                   Specifies that the input is single end reads
    References                      If not specified in the configuration file or you wish to overwrite any of the references.

    Trimming options
    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --clusterOptions              Extra SLURM options, used in conjunction with Uppmax.config
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}
	
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
params.project = false
params.email = false
params.plaintext_email = false

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//Validate inputs	

//if (params.mode != "QC" && params.mode != "characterisation" && params.mode != "complete") {
//	exit 1, "Mode not available. Choose any of <QC, characterisation, complete>"
//}	

if (params.librarylayout != "paired" && params.librarylayout != "single") { 
	exit 1, "Library layout not available. Choose any of <single, paired>" 
}   

if (params.qin != 33 && params.qin != 64) {  
	exit 1, "Input quality offset (qin) not available. Choose either 33 (ASCII+33) or 64 (ASCII+64)" 
}   

//--reads2 can be omitted when the library layout is "single" (indeed it specifies single-end
//sequencing)
if (params.mode != "characterisation" && params.librarylayout != "single" && (params.reads2 == "null") ) {
	exit 1, "If dealing with paired-end reads, please set the reads2 parameters\nif dealing with single-end reads, please set the library layout to 'single'"
}

//--reads1 and --reads2 can be omitted (and the default from the config file used instead) 
//only when mode is "characterisation". Obviously, --reads2 should be always omitted when the
//library layout is single.
if (params.mode != "characterisation" && ( (params.librarylayout == "paired" && (params.reads1 == "null" || params.reads2 == "null")) ||			
 							 			   params.librarylayout == "single" && params.reads1 == "null") ) {
	exit 1, "Please set the reads1 and/or reads2 parameters"
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

refFile = file(params.reference)

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
	set val(pairId), file(R1), file(R2) from ReadPairsToQual

	output:
	set val(pairId), file("${pairId}_dedupe*.fq") into totrim, topublishdedupe

	script:
	"""
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	clumpify.sh -Xmx${maxmem} in1=$R1 in2=$R2 out1=${pairId}_dedupe_R1.fq out2=${pairId}_dedupe_R2.fq \
	qin=$params.qin dedupe subs=0 threads=${task.cpus}
	
	"""
}


/*
 *
 * Step 2: BBDUK: trim + filter (run per sample)
 *
 */

process bbduk {
	tag{ "bbduk" }
	
	input:
   	set val(pairId), file(R1_deduped), file(R2_deduped) from totrim
	file(adapters) from Channel.from( file(params.adapters) )
	file(artifacts) from Channel.from( file(params.artifacts) )
	file(phix174ill) from Channel.from( file(params.phix174ill) )

	output:
	set val(pairId), file("${pairId}_trimmed*.fq") into trimmedreads, todecontaminate, topublishtrim

   	script:
	"""	
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')

	#Quality and adapter trim:
	bbduk.sh -Xmx${maxmem} in=$reads1 in2=$reads2 out=${pairId}_trimmed_R1_tmp.fq \
	out2=${pairId}_trimmed_R2_tmp.fq outs=${pairId}_trimmed_singletons_tmp.fq ktrim=r \
	k=$params.kcontaminants mink=$params.mink hdist=$params.hdist qtrim=rl trimq=$params.phred \
	minlength=$params.minlength ref=$adapters qin=$params.qin threads=${task.cpus} tbo tpe 
	
	#Synthetic contaminants trim:
	bbduk.sh -Xmx${maxmem} in=${pairId}_trimmed_R1_tmp.fq in2=${pairId}_trimmed_R2_tmp.fq \
	out=${pairId}_trimmed_R1.fq out2=${pairId}_trimmed_R2.fq k=31 ref=$phix174ill,$artifacts \
	qin=$params.qin threads=${task.cpus} 

	#Synthetic contaminants trim for singleton read file:
	bbduk.sh -Xmx${maxmem} in=${pairId}_trimmed_singletons_tmp.fq out=${pairId}_trimmed_singletons.fq \
	k=31 ref=$phix174ill,$artifacts qin=$params.qin threads=${task.cpus}

	#Removes tmp files. This avoids adding them to the output channels
	rm -rf ${pairId}_trimmed*_tmp.fq 

	"""
}


/*
 *
 * Step 3: FastQC post-filter and -trim (run per sample)
 *
 */

//This comment is a place holder. The execution of this step is delegated to the process
//Quality Assessment, that is reported at the end of the QC section of this script and
//that runs all the QC assessed tasks.

/*
 *
 * Step 4: Decontamination (run per sample)
 *
 */

process decontaminate {
	tag{ "decon" }
	
	publishDir  workingdir, mode: 'move', pattern: "*_clean.fq.gz"
		
	input:
	set val(pairId), file(infile1), file(infile2), file(infile12) from todecontaminate
	file(refForeingGenome) from Channel.from( file(params.refForeingGenome, type: 'dir') )
	
	output:
	file "*_clean.fq.gz"
	file "${pairId}_clean.fq" into decontaminatedreads, toprofiletaxa, toprofilefunctionreads, topublishdecontaminate

	
	script:
	"""
	
	maxmem=\$(echo ${task.memory} | sed 's/ //g' | sed 's/B//g')
	
	#Decontaminate from foreign genomes
	bbwrap.sh  -Xmx${maxmem} mapper=bbmap append=t in1=$infile1,$infile12 in2=$infile2, null \
	outu=${pairId}_clean.fq outm=${pairId}_cont.fq minid=$params.mind \
	maxindel=$params.maxindel bwr=$params.bwr bw=12 minhits=2 qtrim=rl trimq=$params.phred \
	path=$refForeingGenome qin=$params.qin threads=${task.cpus} untrim quickmatch fast
	
	gzip -c ${pairId}_clean.fq > ${pairId}_clean.fq.gz

	"""
}



/**
	Quality Assessment - STEP 3. Quality control after decontamination. The clean FASTQ 
	file producted by the decontamination step is assessed with FastQC.

	An output similar to that generated by STEP 1 is produced
*/

//This comment is a place holder. The execution of this step is delegated to the process
//Quality Assessment, that is reported below and that runs all the QC assessed tasks.


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
//	QUALITY ASSESSMENT AND VISUALISATION
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

//Creates the correct objects for the quality assessment, by merging the files derived from
//trimming and decontamination and the step number, label and step.
if (params.librarylayout == "paired") {
	trimmedreads2qc = Channel.from('4').combine(trimmedreads.flatMap().merge( Channel.from( ['_R1', '_R2'] ) ){ a, b -> [a, b] }).combine(Channel.from('_trimmedreads'))
} 
else {
	trimmedreads2qc = Channel.from('4').combine(trimmedreads.flatMap()).combine( Channel.from( '' ) ).combine(Channel.from('_trimmedreads'))
}
decontaminatedreads2qc = Channel.from('6').combine(decontaminatedreads).combine( Channel.from( '' ) ).combine(Channel.from('_decontaminatedreads'))

//Creates the channel which performs the QC
toQC = rawreads.mix(trimmedreads2qc, decontaminatedreads2qc) 

//Process performing all the Quality Assessment
process qualityAssessment {
	
	publishDir  workingdir, mode: 'move', pattern: "*.{html,txt}"
	  	
	input:
   	set val(step), file(reads), val(label), val(stem) from toQC

	output:
	file ".log.$step$label" into logQC
	file "${params.prefix}*_fastqc.html" 
	file "${params.prefix}*_fastqc_data.txt" 

	when:
	params.mode == "QC" || params.mode == "complete"

   	script:
	"""	
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Quality Control. [Assessment of read quality] at \$sysdate\" > .log.$step$label
	echo \"File being analysed: $reads\" >> .log.$step$label
	echo \" \" >> .log.$step$label
	
	#Logs version of the software and executed command
	version=\$(fastqc --version) 
	CMD=\"fastqc --quiet --noextract --format fastq --outdir=. --threads ${task.cpus} $reads\"
	
	echo \"Using \$version \" >> .log.$step$label
	echo \"Executing command \$CMD \" >> .log.$step$label
	echo \" \" >> .log.$step$label
	
	#Does QC, extracts relevant information, and removes temporary files
	bash fastQC.sh $reads ${params.prefix}${stem}${label} ${task.cpus} $reads
	
	#Logging QC statistics (number of sequences, Pass/warning/fail, basic statistics, duplication level, kmers)
	base=\$(basename $reads)
	bash logQC.sh \$base ${params.prefix}${stem}${label}_fastqc_data.txt .log.$step$label
				
	#Measures and log execution time			
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"Quality assessment on $reads terminated at \$sysdate (\$exectime seconds)\" >> .log.$step$label
	echo \" \" >> .log.$step$label	
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.$step$label
	echo \" \" >> .log.$step$label	
	"""	
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
//  COMMUNITY CHARACTERISATION 
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

//When doing community characterisation from previously QC'd files, this file should be 
//pushed in the corret channels.
//Please note that the name of the file externally QC'ed should be the same of the 
//one geneared by the YAMP, that is prefix_clean.fq, and should include ALL the reads.
if (params.mode == "characterisation") {
	toprofiletaxa =  Channel.from( file("$workingdir/${params.prefix}_clean.fq") )
	toprofilefunctionreads = Channel.from( file("$workingdir/${params.prefix}_clean.fq") )
}


/**
	Community Characterisation - STEP 1. Performs taxonomic binning and estimates the 
	microbial relative abundancies. MetaPhlAn2 and its databases of clade-specific markers
	are used to infers the presence and relative abundance of the organisms (at the specie/ 
	strain level) that are present in the sample and to estimate their relative abundance.

	Two files are outputted: a tab-separated file reporting the species detected and their
	relative abundance, and a BIOM file, that will be used to evaluate alpha (STEP 8) and
	beta diversity.
*/


process profileTaxa {

	publishDir  workingdir, mode: 'copy', pattern: "*.{biom,tsv}"
	
	input:
	file(infile) from toprofiletaxa
	file(mpa_pkl) from Channel.from( file(params.mpa_pkl) )
	file(bowtie2db) from Channel.fromPath( params.bowtie2db, type: 'dir' )

    output:
	file ".log.7" into log7
	file "${params.prefix}.biom" into toalphadiversity
	file "${params.prefix}_metaphlan_bugs_list.tsv" into toprofilefunctionbugs
	file "${params.prefix}_bt2out.txt" into topublishprofiletaxa

	when:
	params.mode == "characterisation" || params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Community Characterisation. STEP 1 [Taxonomic binning and profiling] at \$sysdate\" > .log.7
	echo \" \" >> .log.7
	
	#If a file with the same name is already present, Metaphlan2 will crash
	rm -rf ${params.prefix}_bt2out.txt
	
	#Defines command for estimating abundances
	CMD=\"metaphlan2.py --input_type fastq --tmp_dir=. --biom ${params.prefix}.biom --bowtie2out=${params.prefix}_bt2out.txt --mpa_pkl $mpa_pkl  --bowtie2db $bowtie2db/$params.bowtie2dbfiles --bt2_ps $params.bt2options --nproc ${task.cpus} $infile ${params.prefix}_metaphlan_bugs_list.tsv\"

	#Logs version of the software and executed command 
	#MetaPhlAn prints on stderr
	version=\$(metaphlan2.py --version 2>&1 >/dev/null | grep \"MetaPhlAn\")
	echo \"Using \$version \" >> .log.7
	echo \"Using BowTie2 database in $params.bowtie2db \" >> .log.7
	echo \" \" >> .log.7
	echo \"Executing command: \$CMD \" >> .log.7
	
	echo \" \" >> .log.7

	#Estimates microbial abundances
	exec \$CMD 2>&1 | tee tmp.log

	#Sets the prefix in the biom file
	sed -i 's/Metaphlan2_Analysis/${params.prefix}/g' ${params.prefix}.biom
	sed -i 's/Metaphlan2_Analysis/${params.prefix}/g' ${params.prefix}_metaphlan_bugs_list.tsv

	#Logs some info
	tree=(kingdom phylum class order family genus species)
	for i in {2..7}
	do
		c=\$(sed '1d' ${params.prefix}_metaphlan_bugs_list.tsv | cut -d\"|\" -f \$i | grep -v \"k__\" | cut -f 1  | sort | uniq | sed '/^\\s*\$/d' | wc -l | cut -d\" \" -f 1)
		echo \"\$c \${tree[((\$i-1))]} found\" >> .log.7
	done

	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"\" >> .log.7
	echo \"STEP 1 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> .log.7
	echo \" \" >> .log.7
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.7
	echo \"\" >> .log.7
	"""
}

/**
	Community Characterisation - STEP 2. Evaluates alpha-diversity, that is the 
	mean species diversity the given sample. Please note that the alpha diversity 
	is the only per-sample measure, so it is the only one evaluated by this module. 
 	
	If a newick tree is provided as input (see QIIME documentation for details), a 
	further and more reliable phylogenetic measure is evaluated (i.e., PD_whole_tree).

	One text file listing the alpha-diversity values, evaluated by means of
	multiple measure, is outputted.
*/

process alphaDiversity {

	publishDir  workingdir, mode: 'move', pattern: "*.{tsv}"
	
	input:
	file(infile) from toalphadiversity
	file(treepath) from Channel.from( file(params.treepath) )
	
    output:
	file ".log.8" into log8
	file "${params.prefix}_alpha_diversity.tsv"
	
	when:
	params.mode == "characterisation" || params.mode == "complete"

	script:
	"""
	#Measures execution time
	sysdate=\$(date)
	starttime=\$(date +%s.%N)
	echo \"Performing Community Characterisation. STEP 2 [Evaluating alpha-diversity] at \$sysdate\" > .log.8
	echo \" \" >> .log.8

	#It checks if the profiling was successful, that is if identifies at least three species
	n=\$(grep -o s__ $infile | wc -l  | cut -d\" \" -f 1)
	if (( n > 3 ))
	then
		#Defines command -- if the tree path is not specified, not all the alpha 
		#measures can be evaluated (that is, PD_whole_tree is skipped)
		if [ $params.treepath == null ]
		then
			CMD=\"alpha_diversity.py -i $infile -o ${params.prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong\"
		else
			CMD=\"alpha_diversity.py -i $infile -o ${params.prefix}_alpha_diversity.tsv -m ace,berger_parker_d,brillouin_d,chao1,chao1_ci,dominance,doubles,enspie,equitability,esty_ci,fisher_alpha,gini_index,goods_coverage,heip_e,kempton_taylor_q,margalef,mcintosh_d,mcintosh_e,menhinick,michaelis_menten_fit,observed_otus,observed_species,osd,simpson_reciprocal,robbins,shannon,simpson,simpson_e,singles,strong,PD_whole_tree -t $treepath\"
		fi
		
		#Logs version of the software and executed command
		version=\$(alpha_diversity.py --version) 
		echo \"Using \$version \" >> .log.8
		if [ $params.treepath == null ]
		then
			echo \"Newick tree not used, PD_whole_tree skipped\" >> .log.8
		else
			echo \"Using Newick tree in $params.treepath\" >> .log.8
		fi
		echo \" \" >> .log.8
		echo \"Executing command: \$CMD \" >> .log.8
		echo \" \" >> .log.8
		
		#Evaluates alpha diversities, redirect is done here because QIIME gets it as an extra parameter
		exec \$CMD 2>&1 | tee tmp.log
	else
		#Also if the alpha are not evaluated the file should be created in order to be returned
		echo \"Not enough classified species detected (N=\$n). Analysis skipped.\" >> .log.8
		touch ${params.prefix}_alpha_diversity.tsv 
	fi
	
	#Measures and log execution time
	endtime=\$(date +%s.%N)
	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
	sysdate=\$(date)
	echo \"\" >> .log.8
	echo \"STEP 2 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> .log.8
	echo \" \" >> .log.8
	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.8
	echo \"\" >> .log.8
	"""
}


/**
	Community Characterisation - STEP 3. Performs the functional annotation using HUMAnN2.
	HUMAnN2 will bypasses the taxomonic profiling step (since it has already
	been performed) and uses the list of specied detected on step 7.
	While the aligners are forced to be Bowtie2 and DIAMOND, the user can
	select the UniRef database to use (UniRef50, UniRef90).

	It outputs several files, and some of them will be removed if keepCCtmpfile is
	set to false. Namely, it creates:
	- three tab-separated files representing the gene families, and the pathways'
	  coverahe and abundancies
	- a SAM file representing the full alignment from Bowtie2 (saved only if keepCCtmpfile
	  is set to true)
    - two tab-separated file representing the reduced aligned reads from both
	  Bowtie2 and DIAMOND (saved only if keepCCtmpfile is set to true)
	- two FASTA file, representing the unaligned reads from both Bowtie2 and
	  DIAMOND (saved only if keepCCtmpfile is set to true)
	- a log of the execution
*/

process profileFunction {

	publishDir  workingdir, mode: 'copy', pattern: "*.{tsv,log}"
	
	input:
	file(cleanreads) from toprofilefunctionreads
	file(metaphlanbuglist) from toprofilefunctionbugs
	file(chocophlan) from Channel.fromPath( params.chocophlan, type: 'dir' )
	file(uniref) from Channel.fromPath( params.uniref, type: 'dir' )
	
    output:
	file ".log.9" into log9
	file "${params.prefix}_HUMAnN2.log"
	file "${params.prefix}_genefamilies.tsv"
	file "${params.prefix}_pathcoverage.tsv"
	file "${params.prefix}_pathabundance.tsv"
	
	//Those may or may be not kept, according to the value of the keepCCtmpfile parameter
	set ("${params.prefix}_bowtie2_aligned.sam", "${params.prefix}_bowtie2_aligned.tsv", "${params.prefix}_diamond_aligned.tsv", "${params.prefix}_bowtie2_unaligned.fa", "${params.prefix}_diamond_unaligned.fa") into topublishhumann2	

	when:
	params.mode == "characterisation" || params.mode == "complete"

	script:
	"""
	#Measures execution time
 	sysdate=\$(date)
 	starttime=\$(date +%s.%N)
 	echo \"Performing Community Characterisation. STEP 3 [Performing functional annotation] with HUMAnN2 at \$sysdate\" > .log.9
 	echo \" \" >> .log.9

	#Defines HUMAnN2 command taking advantages of the MetaPhlAn2's results
	CMD=\"humann2 --input $cleanreads --output . --output-basename ${params.prefix} --taxonomic-profile $metaphlanbuglist --nucleotide-database $chocophlan --protein-database $uniref --pathways metacyc --threads ${task.cpus} --memory-use minimum\"
	
	#Logs version of the software and executed command
	#HUMAnN2 prints on stderr
	version=\$(humann2 --version 2>&1 >/dev/null | grep \"humann2\") 
	echo \"Using \$version \" >> .log.9
	echo \"Using ChocoPhlAn database in $params.chocophlan \" >> .log.9
	echo \"Using UniRef database in $params.uniref \" >> .log.9
	echo \" \" >> .log.9
	echo \"Executing command: \$CMD > ${params.prefix}_HUMAnN2.log\" >> .log.9
	echo \" \" >> .log.9
	
	#Performs functional annotation, redirect is done here because HUMAnN2 freaks out
	#This is  also reported in the log.
	exec \$CMD 2>&1 | tee ${params.prefix}_HUMAnN2.log 

	#If `|| true` is not add, nextflow stops... WTF 
	grep \"Total species selected from prescreen:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Selected species explain\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Unaligned reads after nucleotide alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Total gene families after translated alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	grep \"Unaligned reads after translated alignment:\" ${params.prefix}_HUMAnN2.log >> .log.9 || true
	echo \"More information on HUMAnN2 run are available in the ${params.prefix}_HUMAnN2.log file\" >> .log.9 

	#Some of temporary files (if they exist) may be moved in the working directory, 
	#according to the keepCCtmpfile parameter. Others (such as the bowties2 indexes), 
	#are always removed. Those that should be moved, but have not been created by 
	#HUMAnN2, are now created by the script (they are needed as output for the channel)
	files=(${params.prefix}_bowtie2_aligned.sam ${params.prefix}_bowtie2_aligned.tsv ${params.prefix}_diamond_aligned.tsv ${params.prefix}_bowtie2_unaligned.fa ${params.prefix}_diamond_unaligned.fa)
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

 	#Measures and log execution time
 	endtime=\$(date +%s.%N)
 	exectime=\$(echo \"\$endtime \$starttime\" | awk '{print \$1-\$2}')
 	sysdate=\$(date)
 	echo \"\" >> .log.9
 	echo \"STEP 3 (Community Characterisation) terminated at \$sysdate (\$exectime seconds)\" >> .log.9
 	echo \" \" >> .log.9
 	echo \"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\" >> .log.9
 	echo \"\" >> .log.9
 	"""
}


/**
	CLEANUP 1. Collapses all the logs resulting from the QC in the main one, 
	and removes them.
*/

process logQC {

	input:
	file(tolog)  from logQC.flatMap().mix(log2, log3, log5).toSortedList( { a, b -> a.name <=> b.name } )

	when:
	params.mode == "QC" || params.mode == "complete"

	script:
	"""
	cat $tolog >> $mylog
	"""
}

/**
	CLEANUP 2. Saves the temporary files generate during QC (if the users requested so)
*/
	
	
process saveQCtmpfile {

	publishDir  workingdir, mode: 'copy'
		
	input:
	file (tmpfile) from topublishdedupe.mix(topublishtrim, topublishdecontaminate).flatMap()

	output:
	file "*.fq.gz"

	when:
	(params.mode == "QC" || params.mode == "complete") && params.keepQCtmpfile
		
	script:
	"""
	gzip --force -c $tmpfile > ${tmpfile}.gz
	"""
}

/**
	CLEANUP 3. Collapses all the logs resulting from the community characterisation steps
	in the main one, and removes them.
*/

process logCC {

	input:
	file(tolog) from log7.mix(log8, log9).flatMap().toSortedList( { a, b -> a.name <=> b.name } )
	
	when:
	params.mode == "characterisation" || params.mode == "complete"
		
	script:
	"""
	cat $tolog >> $mylog
	"""
}

/**
	CLEANUP 4. Saves the temporary files generate during the community characterisation 
	(if the users requested so)
*/
	
	
process saveCCtmpfile {

	publishDir  workingdir, mode: 'copy'
		
	input:
	file (tmpfile) from topublishprofiletaxa.mix(topublishhumann2).flatMap()

	output:
	file "$tmpfile"

	when:
	(params.mode == "characterisation" || params.mode == "complete") && params.keepCCtmpfile
		
	script:
	"""
	echo $tmpfile
	"""
}
