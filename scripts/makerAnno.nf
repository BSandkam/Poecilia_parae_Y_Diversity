#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
			                 MAKER ANNOTATION
================================================================================
 #### Author
 Pedro Almeida
--------------------------------------------------------------------------------
*/

// Pipeline version
version = 0.1

/* DEPENDENCIES
RepeatModeler (RepeatModeler, BuildDatabase)
RepeatMasker (RepeatMasker, queryRepeatDatabase.pl)
BUSCO (busco)
SNAP (snap and all associated scripts)
Augustus (augustus)
MAKER (maker and all associated scripts)
BEDtools (bedtools)
seqtk
*/


/* NOTES:
- General:
* fastA file must have read/write access
* several tools do not accept gzipped fasta files, eg. *.fa.gz (gunzip first)

- RepeatMasker / RepeatModeler
* RepeatModeler (optional) requires a minimum of 40000 bp of sequence
* RepeatMasker errors out with long sequence ids (>=50 characters)
  so it is required that the initial fasta ids can be written in less than 20 characters long
* Apparently RepeatModeler does not seem to use all scaffolds/contigs to build a library,
  it samples them, so you might get slightly different results every time you run,
  and might be worth filtering out smaller contigs (e.g. <5000bp)
  so that it does not sample only small contigs for building the libraries.
* PERL5LIB environmental variable must contain the path to the RepeatMasker Perl library

- Augustus
* the `$AUGUSTUS_CONFIG_PATH` path should be writeable without root permissions
*/


/* BEGIN OF SET UP CONFIGURATION VARIABLES */
// --- General options
params.fas = null            // the genome fasta file
params.est = null            // transcripts/EST fasta file
params.prot = null           // protein fasta file
params.prefix = 'makerAnno'  // prefix for the RepeatModeler library
params.out = "$PWD"
params.cpu = 30
// --- RepeatModeler
params.rm_mod = false
params.seqL = 5000           // drop sequences shorter than this length (only for RepeatModeler)
// --- RepeatMasker
params.rm_taxa = null        // the name of the species/taxa for RepeatMasker
params.rm_engine = 'ncbi'    // name of the search engine (abblast/wublast or ncbi (rmblast version))
params.rm_lib = null         // file with a preprepared library for RepeatMasker
params.rm_opts = null        // other options to pass to RepeatMasker
// --- Augustus/BUSCO
params.line = null           // BUSCO's lineage path
params.spc = null            // Augustus species name


params.help = null
if (params.help) {
	log.info ""
	log.info "USAGE: nextflow run makerAnno.nf [--help] [ARGS]"
	log.info ""
	log.info "  Required arguments:"
	log.info "    --fas         FILE  fastA genome file"
	log.info "    --est         FILE  fastA file with transcripts/EST"
	log.info "    --prot        FILE  fastA file with proteins"
	log.info "    --line        DIR   BUSCO's lineage path"
	log.info "    --spc         STR   Augustus species name"
	log.info ""
	log.info "  RepeatModeler arguments:"
	log.info "    --rm_mod      BOOL  run RepeatModeler to create custom library [${params.mod}]"
	log.info "    --seqL        INT   drop sequences shorter than this length [${params.seqL}]"
	log.info ""
	log.info "  RepeatMasker arguments:"
	log.info "    --rm_taxa     STR   name of species/taxa for RepeatMasker"
	log.info "    --rm_engine   STR   name of the search engine [${params.engine}]"
	log.info "    --rm_lib      FILE  file with a preprepared library for RepeatMasker"
	log.info "    --rm_opts     STR   string in quotes with other arguments to pass to RepeatMasker, e.g. ' -nolow -qq'"
	log.info ""
	log.info "  Other optional arguments:"
	log.info "    --prefix      STR   prefix for output files [${params.prefix}]"
	log.info "    --out         DIR   path to output directory [${params.out}]"
	log.info "    --cpu         INT   number of threads [${params.cpu}]"
	log.info ""
	log.info "  Notes:"
	log.info "    * for --rm_opts to work properly it must start with a space"
	log.info "    * fastA file must have short sequence ids (<= 20 characters is recommended)"
	log.info ""
	exit 0
}
/* END OF CONFIGURATION VARIABLES */


if( !params.fas ) exit 1, "No genome fastA file specified!"
if( !params.est ) exit 1, "No EST fastA file specified!"
if( !params.prot ) exit 1, "No file with protein data specified!"
if( !params.line ) exit 1, "No lineage directory specified!"
if( !params.spc ) exit 1, "Missing Augustus species name!"

genome = file(params.fas)
if( !genome.exists() )
	exit 1, "Genome file not found: ${genome}"

lineage = file(params.line)
if( !lineage.exists() )
	exit 1, "Lineage directory not found: ${lineage}"

rm_libfa = null
if (params.rm_lib)
	rm_libfa = file(params.rm_lib)

if( !file(params.est).exists() )
	exit 1, "EST file not found: ${params.est}"

if( !file(params.prot).exists() )
	exit 1, "Protein file not found: ${params.prot}"



/* Header Log info */
log.info "====================================================================="
log.info "             MAKER ANNOTATION pipeline v${version}"
log.info "====================================================================="
log.info "fastA file         : ${params.fas}"
log.info "Prefix             : ${params.prefix}"
log.info "Threads            : ${params.cpu}"
log.info "---------------------------------------------------------------------"
log.info "RepeatModeler"
log.info " run RepeatModeler : ${params.rm_mod}"
log.info "    min seq length : ${params.seqL}"
log.info "---------------------------------------------------------------------"
log.info "RepeatMasker"
log.info "     search engine : ${params.rm_engine}"
log.info "           library : ${params.rm_lib}"
log.info "              taxa : ${params.rm_taxa}"
log.info "  other RM options : ${params.rm_opts}"
log.info "---------------------------------------------------------------------"
log.info "MAKER"
log.info "               est : ${params.est}"
log.info "          proteins : ${params.prot}"
log.info "---------------------------------------------------------------------"
log.info "Augustus/BUSCO"
log.info "           lineage : ${params.line}"
log.info "           species : ${params.spc}"
log.info "---------------------------------------------------------------------"
log.info "Work dir           : ${workDir}"
log.info "Output dir         : ${params.out}"
log.info "====================================================================="
log.info "\n"


/*
------------- FUNCTIONS -------------
*/

/* Function to organise output files */
def organise(def name, def pattern, def istrue, def isfalse,
			 def stdout=null, def stderr=null, def cmdsh=null) {
	if (!stdout) {stdout = '_stdout'}
	if (!stderr) {stderr = '_stderr'}
	if (!cmdsh) {cmdsh = '_cmd'}
	if (name.indexOf(pattern) > 0) istrue
	else if (name ==~ ".command.out") name.replace('.command.out', stdout)
	else if (name ==~ ".command.err") name.replace('.command.err', stderr)
	else if (name ==~ ".command.sh") name.replace('.command.sh', cmdsh)

	else isfalse
}


def fileExists(name){
	// check if a filename name exists. If not exit
	if( !name.exists() ) exit 1, "File not found: ${name}"
}

/*
------------- END OF FUNCTIONS -------------
*/


Channel.fromPath(params.est)
		.splitCsv(sep: '\n', strip: true)
		.map{ it -> file(it[0]) }
		.into{ EST_Check; Ch_EST }
// .subscribe { it -> print it }
EST_Check.subscribe { it -> fileExists(it) }
Ch_EST.collect().set{ EST_2Process }


Channel.fromPath(params.prot)
		.splitCsv(sep: '\n', strip: true)
		.map{ it -> file(it[0]) }
		.into{ Prot_Check; Ch_Prot }
// .subscribe { it -> print it }
Prot_Check.subscribe { it -> fileExists(it) }
Ch_Prot.collect().set{ Prot_2Process }


/* Capture versions from programs into file */
process LogVersions {

	publishDir params.out, mode: 'move'

	output:
	file '_versions.txt'

	script:
	"""
	echo '# RepeatModeler' >> _versions.txt
	RepeatModeler -v >> _versions.txt

	echo '# RepeatMasker' >> _versions.txt
	RepeatMasker -v >> _versions.txt

	echo '# MAKER' >> _versions.txt
	maker -version >> _versions.txt

	echo "# SNAP" >> _versions.txt
	grep version <(snap 2>&1) >> _versions.txt

	echo "# BUSCO" >> _versions.txt
	busco --version >> _versions.txt

	echo "# Augustus" >> _versions.txt
	augustus --version 2>> _versions.txt

	echo "# BEDtools" >> _versions.txt
	bedtools --version >> _versions.txt
	"""
}


OutDir_RM = "${params.out}/01.RM"
if (params.rm_mod) {

	process RunModeler {

		publishDir "${OutDir_RM}/00.modeler", mode: 'copy',
		saveAs:
		{ filename -> organise(filename, ".{fas,fa,stk}", "$filename", "$filename") }

		cpus params.cpu

		input:
		file(genome)

		output:
		file("${params.prefix}-families*") into MODELER_RESULTS
		set file('.command.err'), file('.command.out'), file('.command.sh')

		script:
		Integer threads = params.cpu - 1
		fasta_to_RModel = 'genome2modeler.fas'
		"""
		seqtk seq -l 60 -L ${params.seqL} ${genome} \\
		| fasShortHead.py -l 20 - > ${fasta_to_RModel}

		BuildDatabase \\
			-engine ${params.rm_engine} -name ${params.prefix} ${genome}

		RepeatModeler -srand 12345 \\
			-engine ${params.rm_engine} -pa ${threads} -database ${params.prefix}
		"""
	}

	process RunMaskerFromModeler {

		publishDir "${OutDir_RM}/01.rmasked", mode: 'copy',
		saveAs:
		{ filename -> organise(filename, ".{fas}", "$filename", "$filename") }

		cpus params.cpu

		input:
		file("*") from MODELER_RESULTS

		output:
		file(customLib)
		file("combined.RM.lib")
		file("${outprefix}.*")
		file("${outprefix}.out") into MASKER_OUT
		set file('.command.err'), file('.command.out'), file('.command.sh')

		script:
		Integer threads = params.cpu - 1
		outprefix = genome.fileName.toString()

		customLib = 'rmLib.fa'
		makeCustomLib = "touch ${customLib}"
		if (params.rm_taxa) {
			customLib = params.rm_taxa + '.rmLib.fa'
			makeCustomLib = "queryRepeatDatabase.pl -species ${params.rm_taxa} > ${customLib}"
		}
		"""
		${makeCustomLib}

		cat *-families.fa | awk '{print \$1}' > _tmp.families.fa
		cat _tmp.families.fa ${customLib} > combined.RM.lib

		RepeatMasker \\
			-engine ${params.rm_engine} -pa ${threads} \\
			-noisy \\
			-xsmall \\
			-alignments \\
			-dir . \\
			-lib combined.RM.lib \\
			${params.rm_opts} \\
			${genome}
		"""
	}

} else {

	process RunMasker {

		publishDir "${OutDir_RM}/01.rmasked", mode: 'copy',
		saveAs:
		{ filename -> organise(filename, ".{fas}", "$filename", "$filename") }

		cpus params.cpu

		input:
		file(genome)

		output:
		file("${outprefix}.*")
		file("${outprefix}.out") into MASKER_OUT
		set file('.command.err'), file('.command.out'), file('.command.sh')

		script:
		Integer threads = params.cpu - 1
		outprefix = genome.fileName.toString()

		species = ""
		if (params.rm_taxa) {
			species = "-species ${params.rm_taxa}"
		}

		lib = ""
		if (params.rm_lib) {
			lib = "-lib ${rm_libfa}"
		}

		"""
		RepeatMasker \\
			-engine ${params.rm_engine} -pa ${threads} \\
			-noisy \\
			-xsmall \\
			-alignments \\
			-dir . \\
			${params.rm_opts} \\
			${lib} \\
			${species} \\
			${genome}
		"""
	}
}

process Reformat_RM_GFF {
	publishDir "${OutDir_RM}/02.rm_gff", mode: 'copy',
	saveAs:
	{ filename -> organise(filename, ".gff", "$filename", "$filename") }

	input:
	file(rm_in) from MASKER_OUT

	output:
	file(out)
	file(out_complex)
	file(outGFF) into MASKER_GFF
	set file('.command.err'), file('.command.out'), file('.command.sh')

	script:
	out = rm_in.toString() + '.gff'
	out_complex = out - ~/\.gff$/ + '.complex.gff'
	outGFF = out_complex - ~/\.gff$/ + '.reformat.gff'
	"""
	# convert out to gff3
	rmOutToGFF3.pl ${rm_in} > ${out}

	# separate complex from simple repeats. Only the former will be used in MAKER
	# - https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2
	# - https://groups.google.com/forum/#!topic/maker-devel/patU-l_TQUM
	grep -v -e "Satellite" -e ")n" -e "-rich" ${out} > ${out_complex}

	# reformat to work with MAKER
	fixgff3_script.pl ${out_complex} > ${outGFF}
	"""
}

base_rnd1 = "${params.prefix}.round1"
process MAKER_Round1 {
	publishDir "${params.out}/02.MAKER_round1", mode: 'copy',
	saveAs:
	{ filename -> organise(filename, "{.ctl}", "$filename", "$filename") }

	cpus params.cpu

	input:
	file(genome)
	file(est) from EST_2Process
	file(prot) from Prot_2Process
	file(rm_gff) from MASKER_GFF

	output:
	file("*.ctl")
	file("${base_rnd1}.maker.output") into (
										MAKER_RND1_DIR2SNAP, MAKER_RND1_DIR2AUG,
										MAKER_RND1_RND2)
	set file('.command.err'), file('.command.out'), file('.command.sh')

	script:
	est_list = est.join(',')
	prot_list = prot.join(',')
	ctl = "maker_opts_round1.ctl"
	"""
	echo "genome=${genome}" > ${ctl}
	echo "organism_type=eukaryotic" >> ${ctl}
	echo "est=${est_list}" >> ${ctl}
	echo "protein=${prot_list}" >> ${ctl}
	echo "model_org=simple" >> ${ctl}
	echo "rm_gff=${rm_gff}" >> ${ctl}
	echo "est2genome=1" >> ${ctl}
	echo "protein2genome=1" >> ${ctl}
	echo "cpus=${task.cpus}" >> ${ctl}

	maker -BOPTS
	maker -EXE

	maker -base ${base_rnd1} ${ctl} maker_bopts.ctl maker_exe.ctl

	cd ${base_rnd1}.maker.output
		gff3_merge -s -d ${base_rnd1}_master_datastore_index.log > ${base_rnd1}.all.maker.gff
		fasta_merge -d ${base_rnd1}_master_datastore_index.log
		# GFF w/o the sequences
		gff3_merge -n -s -d ${base_rnd1}_master_datastore_index.log > ${base_rnd1}.all.maker.noseq.gff
	cd -
	"""
}


OutDir_Training_round1 = "${params.out}/03.Training_Round1"
process Train_SNAP_Round1 {
	publishDir "${OutDir_Training_round1}/01.SNAP", mode: 'copy',
	saveAs:
	{ filename -> organise(filename, ".{hmm,log}", "$filename", "$filename") }

	input:
	file(maker1dir) from MAKER_RND1_DIR2SNAP

	output:
	file("${params.prefix}.snap.rnd1.hmm") into SNAP_RND1_HMM
	file("*.log")
	set file('.command.err'), file('.command.out'), file('.command.sh')

	script:
	"""
	# export 'confident' gene models from MAKER
	maker2zff -o 0 -a 0 -e 0 -c 0 -x 0.5 -l 50 -d ${maker1dir}/${base_rnd1}_master_datastore_index.log

	fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
	fathom genome.ann genome.dna -validate > validate.log 2>&1

	# collect the training sequences and annotations, plus 1000 surrounding bp for training
	fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
	fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

	mkdir params && cd params
	forge ../export.ann ../export.dna > ../forge.log 2>&1
	cd ..

	hmm-assembler.pl ${params.prefix}.snap.rnd1 params > ${params.prefix}.snap.rnd1.hmm
	"""
}


process Train_Augustus_Round1 {
	publishDir "${OutDir_Training_round1}/02.Augustus", mode: 'copy',
	saveAs:
	{ filename -> organise(filename, ".fa", "$filename", "$filename") }
	cpus params.cpu

	input:
	file(genome)
	file(maker1dir) from MAKER_RND1_DIR2AUG

	output:
	file("run_${busco_out}/augustus_output/retraining_parameters") into AUGUSTUS_RND1_PARS
	file(transcripts)
	set file('.command.err'), file('.command.out'), file('.command.sh')

	script:
	transcripts = "${base_rnd1}.all.maker.transcripts1000.fa"
	busco_out = "${params.prefix}.busco.rnd1"
	"""
	cat ${maker1dir}/${base_rnd1}.all.maker.noseq.gff \\
	| awk -v OFS="\\t" '{ if (\$3 == "mRNA") print \$1,\$4,\$5 }' \\
	| awk -v OFS="\\t" \\
		'{ if (\$2 < 1000) print \$1,"0",\$3+1000; else print \$1,\$2-1000,\$3+1000 }' \\
	| bedtools getfasta -fi ${genome} -bed - -fo ${transcripts}

	busco -i ${transcripts} -o ${busco_out} \\
		-m genome -c ${task.cpus} --long \\
		--augustus_parameters='--progress=true' \\
		-l ${params.line} -sp ${params.spc}
	"""
}

/*
The following process simply renames Augustus parameter files
and copies them into $AUGUSTUS_CONFIG_PATH.
Note that if a folder named ${params.prefix}.rnd1 already exists in $AUGUSTUS_CONFIG_PATH
it will be removed and replaced by a new one.
*/
augSpName1 = "${params.prefix}_rnd1"
process Augustus_Round1_Process {

	input:
	file(augustusDir) from AUGUSTUS_RND1_PARS

	output:
	file(checkpoint) into AUGUSTUS_RND1_CHECKPOINT

	script:
	checkpoint = '_checkpoint.augustus1'
	"""
	cp -p ${augustusDir}/* .

	_prefix=\$(ls *_parameters.cfg | sed 's/_parameters.cfg//')

	for f in \${_prefix}_*; do
		sed -i '' "s/\${_prefix}/${augSpName1}/g" \$f
		mv \$f \${f/\${_prefix}/${augSpName1}}
	done

	rm -rf \$AUGUSTUS_CONFIG_PATH/species/${augSpName1}
	mkdir \$AUGUSTUS_CONFIG_PATH/species/${augSpName1}
	cp ${augSpName1}_* \$AUGUSTUS_CONFIG_PATH/species/${augSpName1}

	touch ${checkpoint}
	"""
}


base_rnd2 = "${params.prefix}.round2"
process MAKER_Round2 {
	publishDir "${params.out}/04.MAKER_round2", mode: 'copy',
	saveAs:
	{ filename -> organise(filename, "{.ctl}", "$filename", "$filename") }

	cpus params.cpu

	input:
	file(genome)
	file(maker1dir) from MAKER_RND1_RND2
	file(snap1_hmm) from SNAP_RND1_HMM
	file(_) from AUGUSTUS_RND1_CHECKPOINT

	output:
	file("*.ctl")
	file("${base_rnd2}.maker.output") into (
										MAKER_RND2_DIR2SNAP, MAKER_RND2_DIR2AUG,
										MAKER_RND2_RND3)
	set file('.command.err'), file('.command.out'), file('.command.sh')

	script:
	ctl = "maker_opts_round2.ctl"
	"""
	# transcript alignments from 1st round
	cat ${maker1dir}/${base_rnd1}.all.maker.noseq.gff \\
	| awk '{ if (\$2 == "est2genome") print \$0 }' > ${base_rnd1}.all.maker.est2genome.gff
	# protein alignments from 1st round
	cat ${maker1dir}/${base_rnd1}.all.maker.noseq.gff \\
	| awk '{ if (\$2 == "protein2genome") print \$0 }' > ${base_rnd1}.all.maker.protein2genome.gff
	# repeat alignments from 1st round
	cat ${maker1dir}/${base_rnd1}.all.maker.noseq.gff \\
	| awk '{ if (\$2 ~ "repeat") print \$0 }' > ${base_rnd1}.all.maker.repeats.gff


	echo "genome=${genome}" > ${ctl}
	echo "organism_type=eukaryotic" >> ${ctl}
	echo "est=" >> ${ctl}
	echo "est_gff=${base_rnd1}.all.maker.est2genome.gff" >> ${ctl}
	echo "protein=" >> ${ctl}
	echo "protein_gff=${base_rnd1}.all.maker.protein2genome.gff" >> ${ctl}
	echo "model_org=" >> ${ctl}
	echo "rm_gff=${base_rnd1}.all.maker.repeats.gff" >> ${ctl}
	echo "snaphmm=${snap1_hmm}" >> ${ctl}
	echo "augustus_species=${augSpName1}" >> ${ctl}
	echo "est2genome=0" >> ${ctl}
	echo "protein2genome=0" >> ${ctl}
	echo "split_hit=20000" >> ${ctl}
	echo "cpus=${task.cpus}" >> ${ctl}

	maker -BOPTS
	maker -EXE

	maker -base ${base_rnd2} ${ctl} maker_bopts.ctl maker_exe.ctl

	cd ${base_rnd2}.maker.output
		gff3_merge -s -d ${base_rnd2}_master_datastore_index.log > ${base_rnd2}.all.maker.gff
		fasta_merge -d ${base_rnd2}_master_datastore_index.log
		# GFF w/o the sequences
		gff3_merge -n -s -d ${base_rnd2}_master_datastore_index.log > ${base_rnd2}.all.maker.noseq.gff
	cd -
	"""
}


workflow.onComplete {
	log.info "\n"
	log.info "========================================="
	log.info "           MAKER ANNOTATION              "
	log.info "========================================="
	log.info "Session ID     : ${workflow.sessionId}"
	log.info "Res. previous  : ${workflow.resume}"
	log.info "Completed at   : ${workflow.complete}"
	log.info "Duration       : ${workflow.duration}"
	log.info "Success        : ${workflow.success}"
	log.info "Workdir        : ${workflow.workDir}"
	log.info "Outdir         : ${params.out}"
	log.info "Exit status    : ${workflow.exitStatus}"
	log.info "Error report   : ${workflow.errorReport ?: '-'}"
	log.info "========================================="
	log.info "\n"
}
