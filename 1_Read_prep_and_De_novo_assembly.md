**Goals:**
  1) To generate de novo genomes using the 10X supernova pipeline
  2) To prepare the 10X raw sequencing reads for downstream alongside straight Illumina Reads**
  3) To prepare the straight Illumina reads for downstream analysis


# **Step 1) Generate de novo genomes using the Supernova (10X Genomics) pipeline**

## 1a) Run Supernova with reads that have gone through the 10X Chromium linked-read pipeline

Parameters:
  supernova run (v2.1.1)
    --id = ID for the run that is used to name the output
    --outprefix = prefix to use for fasta files created
    --maxreads==all = uses all sequencing reads available. Can be downsampled if the number of reads is making the run difficult. Downsampling did not improve the de novo genomes created using P. parae sequence data, therefore all reads were used for each individual. Uses a large amount of RAM - each run required ~1.3 TB RAM to run.
    --accept-extreme-coverage = flag to override the early exit when using large number of reads
    --fastqs = path to directory that contains the 10X linked-read directory structure with the fastq files
    --sample = prefix of the sequence files for the sample that are located in the fastqs directory indicated above

```bash
supernova run --id=Sample01 --maxreads=all --accept-extreme-coverage --fastqs=../[reads_location] --sample=[reads_ID]
```


## 1b) Generate Pseudohap2 style output from 10X Supernova pipeline

Parameters:
  supernova mkoutput (v2.1.1)
    --asmdir = path to directory with the final outputs from the supernova run
    --outprefix = prefix to use for fasta files created
    --style=pseudohap2 = fasta output that arbitrarily chooses sequence differences that may be coming from two chromosomes (thus is unphased) to make a single fasta sequence (akin to generating a single fasta sequence for a genome). See description and illustrations here: https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/output/generating

```bash
supernova mkoutput --asmdir=../Sample01/outs/assembly --outprefix=Sample01_pseudohap --style=pseudohap2 --headers=short --index
```


## 1c) Generate Megabubble style output from 10X Supernova pipeline

Parameters:
  supernova mkoutput (v2.1.1)
    --asmdir = path to directory with the final outputs from the supernova run
    --outprefix = prefix to use for fasta files created
    --style=megabubbles = fasta output that arbitrarily flattens very minor unresolved sequence differences (often technical error) but does not attempt to flatten major sequence differences that are likely coming from two chromosomes. See descriptions and illustrations here: https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/output/generating

```bash
supernova mkoutput --asmdir=../Sample01/outs/assembly --outprefix=Sample01_megabubble --style=megabubbles --headers=short
```


# **Step 2) Prepare 10X raw sequencing reads for downstream alongside straight Illumina Reads**

## 2a) Run Longranger basic (10X Genomics) pipeline to strip the barcodes from the linked-read pipeline from the Illumina reads

Parameters:
  longranger (v2.2.2)
    --sample = sample prefix from sequencer
    --id = sample name
    --fastqs = path to directory that contains the 10X linked-read directory structure with the fastq files

```bash
longranger basic --sample=[sampleID_from_sequencer] --id=Sample01 --fastqs=../[reads_location]
```


## 2b) Run Trimmomatic on reads from 10X linked-read pipeline

Parameters:
  trimmomatic (v0.36)
    PE = Indicates paired end reads
    -threads = 25 threads
    -phred33 = indicates phred33 quality score
    ILLUMINACLIP = indicates the Illumina library to be clipped from reads
    LEADING:3 = trims reads from start until base of quality 3
    TRAILING:3 = trims reads back from end until base of quality 3
    SLIDINGWINDOW:4:15 = slides window along read to trim if internal quality is too low
    MINLEN:50 = minimum length of read required after trimming (otherwise read is discarded)
*This step uses the script deinterleave_fastq.sh available at https://gist.github.com/3521724 to de-interleave the fastq files (from supernova basic) so that reads can be treated the same whether they went through the 10X linked-read pipeline or straight Illumina sequencing*

```bash
gzip -dc ../from_longranger/Sample01/outs/barcoded.fastq.gz | deinterleave_fastq.sh Sample01_for.fastq Sample01_rev.fastq
java -jar trimmomatic-0.36.jar PE -threads 12 -phred33 Sample01_for.fastq Sample01_rev.fastq Sample01_for_paired.fastq.gz Sample01_for_unpaired.fastq.gz Sample01_rev_paired.fastq.gz Sample01_rev_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```


# **Step 3) Prepare straight Illumina reads for downstream analysis**

Parameters:
  trimmomatic (v0.36)
    PE = Indicates paired end reads
    -threads = 25 threads
    -phred33 = indicates phred33 quality score
    ILLUMINACLIP = indicates the Illumina library to be clipped from reads
    LEADING:3 = trims reads from start until base of quality 3
    TRAILING:3 = trims reads back from end until base of quality 3
    SLIDINGWINDOW:4:15 = slides window along read to trim if internal quality is too low
    MINLEN:50 = minimum length of read required after trimming (otherwise read is discarded)

```bash
java -jar trimmomatic-0.36.jar PE -threads 25 -phred33 [forward_reads.fastq.gz] [reverse_reads.fastq.gz] [for_paired.fastq.gz] [for_unpaired.fastq.gz] [rev_paired.fastq.gz] [rev_unpaired.fastq.gz] ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```
