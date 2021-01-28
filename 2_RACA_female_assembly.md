**Goal:**
  Use the Reference Assisted Chromatin Assembly (RACA) pipeline to order and orient *de novo* female scaffolds by chromosome

RACA needs a few datasets to run:
- A *de novo* genome done to scaffold level of Poecilia parae (target species)
- Paired End Illumina sequences with short inserts (~500bp) from additional individuals (the more the better - we used 9 additional female samples)
- Mate Pair Illumina sequences with large inserts (several kb+) from additional individuals. Rather than generating additional Mate Pair libraries we generated "pseudo-mate pair libraries" (with 2kb and 15kb inserts) from the *de novo* genomes of the additional female samples sequenced with the 10X linked-reads and assembled with Supernova (10X Genomics) described previously.
- Chromosome level genome of a closely related species: Here we use Xiphophorus helleri
- Chromosome level genome of an outgroup species (older than split between above species' genome and your target): Here we use Oryzias latipes
- A phylogeny that contains just the target, the reference, and the outgroup (with node ages)


Software used:
- [KentUtils](https://github.com/ENCODE-DCC/kentUtils)
- [LastZ](http://www.bx.psu.edu/~rsharris/lastz/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [RACA](http://bioinfo.konkuk.ac.kr/RACA/) I used RACA v1.1 which is on this site here. I think it's got a bit more bug fixes than v1.0 (which I think is what is on GitHub?)


Scripts used:
- sam2raca_convert.pl
- size_selecting.pl
- mate pair maker


# Steps in the RACA pipeline:

1. Prepare RACA input files
   * a. Prepare best *de novo* female *Poecilia parae* genome to run RACA on (_aka 'target'_)
   * b. Map the Illumina paired end reads (short inserts) and the Illumina mate pair reads (long inserts) to your Parae female genome
   * c. Generate RACA-ready 'reads' files from your mapped reads

2. Generate alignments of Xiphophorus helleri (a close relative) reference genome to (1) Parae female (target) genome assembly and (2) Oryzias latipes (an outgroup ancestral species) reference genome
   * 2a. Generate assembly sizes for (1) scaffolds of Parae, (2) scaffolds of close ancestor, and (3) each of the outgroup relative's chromosomes.
   * 2b. Convert to .2bit
   * 2c. Run LastZ between each of the close relative's chromosomes
   * 2d. Convert lastz fasta to .psl
   * 2e. Produce Parae-Reference Chains
   * 2f. Produce Parae-Reference Nets
   * 2g. Produce Outgroup-Reference Chains
   * 2h. Produce Outgroup-Reference Nets
   * 2i. Move Chain and Net files into file structure for RACA

3. Adjusting parameters and running RACA
   * 3a. Update and save a 'config.SFs'
   * 3b. Update and save a 'insertsize_sd.txt'
   * 3c. Update and save a 'params.txt'
   * 3d. Update and save a 'readmapping_lib.txt'
   * 3e. Save a 'reliable_adjs.txt'
   * 3f. Update and save a 'tree.txt'
   * 3g. Run RACA




# **Step 1) Prepare RACA input files**

## 1a - Prepare your *de novo* genome of the species you want to run RACA for
### 1a.1) Make reference file of just the 10kb plus scaffolds.
I used only scaffolds that were at least 10kb long because smaller scaffolds give RACA some trouble as they tend to be highly repetitive and polymorphic sequence which creates ambiguities and reduces RACA accuracy.

*Note: This requires the following script called size_selecting.pl*
```Perl
## size_selecting.pl
# This is a script to size select a fasta file.
# Invoke as >   perl size_selecting.pl {min_length} {fasta_file} > {new_file_name.fasta}
#!/usr/bin/perl
use strict;
use warnings;

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}
```

Parameters
  size_selecting.pl
    10000 = only keep scaffolds 10kb+

```bash
perl size_selecting.pl 10000 Parae_female.fasta > Parae_female_10kb.fasta
```

### 1a.2) Take the output file from the Supernova run and make all the names 'scaffold' followed by unique number.
This will depend upon the format of your *de novo* genome file. RACA needs each scaffold to be called '>scaffold#' (eg. '>scaffold12'). Output from supernova resulted in scaffolds that were just called by a number. Therefore I had to change the name header of each scaffold to include the word 'scaffold'.
*Note: It's really important you get your names fixed here otherwise the rest of the pipeline will break at the end.*

```bash
sed 's/>/>scaffold/g' [Parae_female_10kb.fasta] > Parae_female_10kb_scaffolds.fasta
```

### 1a.3) Generate index for the 10-kb filtered genome
You'll need an index for your *de novo* genome so you can align things to it.

Parameters
  bowtie2-build (v.2.2.9)
    -f = file of scaffolds to index

- [ ] Index your genome.
```bash
bowtie2-build -f Parae_female_10kb_scaffolds.fasta Parae_female_10kb_scaffolds.index
```


## 1b) Map reads to assembly
Map both Paired-end reads with short inserts and mate pair reads with long inserts to Parae female scaffolds.

### 1b.1) Map paired-end short insert reads for each of additional Parae female samples to Parae female genome.
Parameters
  bowtie2 (v.2.2.9)
    -q = input files are in fastq format
    -p 12 = run with 12 threads
    --no-discordant = supresses discordant alignments for paired reads
    -X 800 = maximum fragment insert expected
    --fr = expected read orientation is forward-reverse (for paired-end reads)
    -x = index filename
    -1 = file of forward reads
    -2 = file of reverse reads
    -S = name of sam file to save output to

```bash
bowtie2 -q -p 12 --no-discordant -X 800 --fr -x Parae_female_10kb_scaffolds.index -1 sample1_forward_paired.fastq -2 sample1_forward_paired.fastq -S Parae_female_illum1.sam
```

### 1b.2a) Generate 2kb and 15kb pseudo mate-pair long-insert reads from the de novo assembled linked-reads of P. parae.

*Note: this use the python script mate_pair_maker.py to create pseudo mate-pairs from de novo scaffolds - that is to say that it takes scaffolds that are longer than the requested insert size and generates a fastq file of forward and reverse reads with the designated insert size and in the reverse forward orientation expected by mate pair library data*

Parameters
  -f = input fasta file of de novo assembly
  -in = insert size
  -seq = the length of forward and reverse sequences generated
  -step = the step size between reads generated
  -out = the outfile basename

```bash
#For each female sample with de novo sequence
#Generate 2kb inserts
python3.5 mate_pair_maker.py -f [female1_denovo_assembly] -in 2000 -seq 150 -step 12 -out [sample1_2kb_pseudo-MatePairs]
#Generate 2kb inserts
python3.5 mate_pair_maker.py -f [female1_denovo_assembly] -in 15000 -seq 150 -step 12 -out [sample1_2kb_pseudo-MatePairs]
```


### 1b.2b) Map the 2kb and 15kb pseudo mate-pair long-insert reads for each additional sample to Parae female genome.

Parameters
  bowtie2 (v.2.2.9)
    -q = input files are in fastq format
    -p 12 = run with 12 threads
    --no-discordant = supresses discordant alignments for paired reads
    -I 1200 = the minimum insert span size
    -X 4000 = maximum insert size expected
    --rf = expected read orientation is reverse-forward (for Mate-pair library reads)
    -x = index filename
    -1 = file of forward reads
    -2 = file of reverse reads
    -S = name of sam file to save output to
*Note: the -I option is the minimum insert span, and the -X is the maximum insert span. If you have a true Mate-Pair library with a range you can estimate this as I = (2x read_length) + (0.5x insert_length) and X = I + 1000*
```bash
bowtie2 q -p 12 --no-discordant -I 1200 -X 4000 --rf -x genome_assembly -1 sample1_2kb_pseudo-MatePairs_forward.fastq -2 sample1_2kb_pseudo-MatePairs_reverse.fastq -S Parae_female_sample1_2kb_mp_bowtie.sam

bowtie2 q -p 12 --no-discordant -I 12000 -X 20000 --rf -x genome_assembly -1 sample1_15kb_pseudo-MatePairs_forward.fastq -2 sample1_15kb_pseudo-MatePairs_reverse.fastq -S Parae_female_sample1_15kb_mp_bowtie.sam
```


## 1c) Convert the sam files to a RACA ready format
The sam2raca_convert.pl script is part of the RACA package, it takes a sam file and generates a mapping file, and a mapping stats file. The 3rd option you need to specify is whether the sam is (A) an 'innie' - this is when the direction of paired end reads face one another (as in normal Illumina paired end sequencing ->  <- ), or if the reads are (B) 'outie' - this is when the direction of paired end reads face away from one another (as in normal Mate Pair library generation <- -> ).

*Note: Run the script for each sample and each sequencing style (paired end and mate pair sequencing*
```bash
sam2raca_convert.pl Parae_female_illum1.sam Parae_female_illum1_mapping_file Parae_female_illum1_mapping_stats innie

mkdir ../Running_RACA
mkdir ../Running_RACA/reads
for file in *_mapping_file ; do mv $file ../Running_RACA/reads ; done
```

# **Step 2) Generate alignments of Parae female genome assembly and the Xiphophorus helleri (close relative)**
*Note: The next several steps utilize scripts in the [KentUtils](https://github.com/ENCODE-DCC/kentUtils) software suite*

### 2a. Generate assembly sizes for (1) scaffolds of Poecilia parae female (target), (2) chromosomes of Xiphophorus helleri (close relative), and (3) chromosomes of Oryzias latipes (outgroup relative).

```bash
faSize Parae_female_10kb_scaffolds.fasta -detailed > Parae_female_10kb_scaffolds.assembly.size
cp Parae_female_10kb_scaffolds.assembly.size ../Running_RACA/
cp Parae_female_10kb_scaffolds.fasta ../Running_RACA/
faSize Xiphophorus.fa -detailed > Xiphophorus.assembly.size
```

*For the reference genome of the close relative each chromosome must be it's own file and the fasta header for each chromosome must start with "chr#"" (eg. chr1, chr2, etc).*

```bash
for file in *.fa ; do faSize $file -detailed > $file.assembly.size ; done
```

### 2b. Use the faToTwoBit in the kentUtils package to convert the fasta files to 2bit format

```bash
#Poecilia parae (sample P32)
faToTwoBit Parae_female_10kb_scaffolds.fasta Parae_female_assembly.2bit
#Oryzias latipes
faToTwoBit olat.fa olat.2bit
# Run in folder with individual files of close relative's chromosomes called 'xiph_LG1.fa' etc. Xiphophorus helleri has 24 chromosomes so that is how it is written and this is used throughout- should be adjusted accordingly.
for i in {1..24} ; do faToTwoBit xiph_LG$i.fa xiph_LG$i.2bit ; done
```

### 2c. Run LastZ to obtain pairwise alignments between each reference chromosome of Xiphophorus helleri (close relative) and Poecilia parae (target), and Oryzias latipes (outgroup relative) genomes.

Parameters
  lastz (v.1.04.00)
    C=0 = same as --nochain --gapped
    H=2000 = inner score 2000
    M=50 = mask any position in target occuring more than 50 times
    ambiguous=iupac = treat any ambiguous IUPAC-IUB character as a completely ambiguous nucleotide
    --markend = marks last line of file to indicate lastz ran to completion
    --rdotplot = save dotplot of alignment as filename

```bash
#Run for outgroup (Oryzias latipes)
for i in {1..24} ; do lastz xiph_LG$i.2bit olat.2bit C=0 H=2000 M=50 --ambiguous=iupac --output=lastZ_xiph_olat_LG$i.fa --markend --rdotplot=lastZ_xiph_olat_LG$i ; done

#Run for target (P. parae)
for i in {1..24} ; do lastz xiph_LG$i.2bit Parae_female_assembly.2bit C=0 E=30 H=2000 K=3000 L=3000 O=400 M=50 --ambiguous=iupac --output=lastz_xiph_Parae_female_LG$i.fa --markend --rdotplot=lastz_xiph_Parae_female_LG$i ; done
```


* The resulting LastZ alignment fasta files are then run through the [UCSC Chains and Nets](https://www.pnas.org/content/100/20/11484.short) pipeline. The UCSC chain and net nucleotide alignment representation formats can accommodate the existence of sequence duplications, deletions, transpositions and inversions what makes them a helpful tool in the study of genome evolution. “Chains” represent an ordered sequence of traditional pairwise nucleotide alignments separated by large gaps that can be present in both species genomes. “Nets” represent the chains that cover a same sequence region in a hierarchical manner, starting with the longest and highest synteny-scoring chain for that genome region. RACA will use these files to define SFs, merging co-linear alignments. More details on the Chains and Nets can be found here: http://genomewiki.ucsc.edu/index.php/Chains_Nets.

*The following steps are all with commands from KentUtils*

## 2d. Convert lastz fasta to .psl format

```bash
for i in {1..23} ; do lavToPsl lastz_xiph_Parae_female_LG$i.fa lastz_xiph_Parae_female_LG$i.psl ; done

# Move files to respective LG folders
for i in {1..24} ; do mkdir LG$i ; done
for i in {1..24} ; do mv lastz_xiph_Parae_female_LG$i.psl LG$i\/ ; done
```

### 2e. Produce Parae-Xiphophorus (Target-Reference) Chains

Parameters (all part of KentUtils package)
  axtChain
    -minScore=1000 = minimum score for a chain
    -linearGap=medium = specifies the linear gap to use (medium ~ mouse/human linear gap costs)
    -verbose=0 = run quietly
    -psl = use .psl format

```bash
# Run axtChain to chain together alignments
for i in {1..24} ; do axtChain -minScore=1000 -linearGap=medium -verbose=0 -psl LG$i\/lastz_Parae_female_xiph_LG$i.psl /xiph_LG$i.2bit Parae_female_assembly.2bit LG$i/Parae_female_xiph_LG$i.chain ; done

# Run chainAntiRepeat to get rid of chains that are primarily the results of repeats and degenerate DNA
for i in {1..24} ; do chainAntiRepeat xiph_LG$i.2bit Parae_female_assembly.2bit LG$i/Parae_female_xiph_LG$i.chain LG$i/Parae_female_xiph_antirepeat_LG$i.chain ; done

# Run chainSort to sort chains by score
for i in {1..24} ; do chainSort LG$i/Parae_female_xiph_antirepeat_LG$i.chain LG$i/Parae_female_xiph_sort_LG$i.chain ; done
```


### 2f. Produce Parae-Xiphophorus (Target-Reference) Nets

```bash
# Run chainPreNet to remove chains that don't have a chance of being netted
for i in {1..24} ; do chainPreNet LG$i/Parae_female_xiph_sort_LG$i.chain xiph_chromosome$i.fa.assembly.size Parae_female_assembly.size LG$i/Parae_female_xiph_prenet_LG$i.chain ; done

# Run chainNet to make alignment nets out of chains
for i in {1..24} ; do chainNet LG$i/Parae_female_xiph_prenet_LG$i.chain xiph_chromosome$i.fa.assembly.size Parae_female_10kb_scaffolds.assembly.size LG$i/xiph_Parae_female_LG$i.net LG$i/Parae_female_xiph_LG$i.net ; done

# Run netSyntenic to add synteny info to nets
for i in {1..24} ; do netSyntenic LG$i/xiph_Parae_female_LG$i.net LG$i/xiph_Parae_female_noClass_LG$i.net ; done
```


*Note: Now do essentially the same thing for the outgroup and the closely related reference*
### 2g. Produce Outgroup-Reference Chains

```bash
# Convert to psl format
for i in {1..24} ; do lavToPsl lastZ_xiph_olat_LG$i.fa lastZ_xiph_olat_LG$i.psl ; done

# Move files to respective LG folders
for i in {1..24} do mkdir LG$i ; done
for i in {1..24} ; do mv lastZ_xiph_olat_LG$i.psl LG$i\/ ; done

# Run axtChain to chain together alignments
for i in {1..24} ; do axtChain -minScore=1000 -linearGap=medium -verbose=0 -psl LG$i\/lastZ_xiph_olat_LG$i.psl xiph_LG$i.2bit olat.2bit LG$i/olat_xiph_LG$i.chain ; done

# Run chainAntiRepeat to get rid of chains that are primarily the results of repeats and degenerate DNA
for i in {1..24} ; do chainAntiRepeat xiph_LG$i.2bit olat.2bit LG$i/olat_xiph_LG$i.chain LG$i/olat_xiph_antirepeat_LG$i.chain ; done

# Run chainSort to sort chains by score
for i in {1..24} ; do chainSort LG$i/olat_xiph_antirepeat_LG$i.chain LG$i/olat_xiph_sort_antirepeat_LG$i.chain ; done
```

### 2h. Produce Outgroup-Reference Nets

```bash
# Run chainPreNet to remove chains that don't have a chance of being netted
for i in {1..24} ; do chainPreNet LG$i/olat_xiph_sort_antirepeat_LG$i.chain xiph_chromosome$i.fa.assembly.size olat.assembly.size LG$i/olat_xiph_prenet_LG$i.chain ; done

# Run chainNet to make alignment nets out of chains
for i in {1..24} ; do chainNet LG$i/olat_xiph_prenet_LG$i.chain xiph_chromosome$i.fa.assembly.size olat.assembly.size LG$i/xiph_olat_LG$i.net LG$i/olat_xiph_LG$i.net ; done

# Run netSyntenic to add synteny info to nets
for i in {1..24} ; do netSyntenic LG$i/xiph_olat_LG$i.net LG$i/xiph_olat_noClass_LG$i.net ; done
```

### 2i. Move chain and net files into file structure for RACA
*Note: RACA requires a very specific directory structure*

```bash
mkdir ../Running_RACA/xiph.chain.net
mkdir ../Running_RACA/xiph.chain.net/xiph
mkdir ../Running_RACA/xiph.chain.net/xiph/ppar
mkdir ../Running_RACA/xiph.chain.net/xiph/ppar/chain
mkdir ../Running_RACA/xiph.chain.net/xiph/ppar/net
for i in {1..24} ; do cp LG$i/parae_xiph_sort_LG$i.chain Running_RACA/xiph.chain.net/xiph/ppar/chain/chr$i.chain ; done
for i in {1..24} ; do cp LG$i/Parae_female_xiph_noClass_LG$i.net ../Running_RACA/xiph.chain.net/xiph/ppar/net/chr$i.net ; done
for i in {1..24} ; do cp LG$i/olat_xiph_sort_antirepeat_LG$i.chain ../Running_RACA/xiph.chain.net/xiph/olat/chain/chr$i.chain ; done
for i in {1..24} ; do cp LG$i/xiph_olat_noClass_LG$i.net ../Running_RACA/xiph.chain.net/xiph/olat/net/chr$i.net ; done
```


# **Step 3) Run RACA**
*Note: there are several configuration files that need to be updated and saved before running. Adjust accordingly*
### 3a. Update the 'config.SFs' file
*Following is the list of things that need to be adjusted in the config file for a specific run*
- [ ] Make the 'netdir' to be the full path to /Running_RACA/xiph.chain.net
- [ ] Make the 'chaindir' to be the full path to /Running_RACA/xiph.chain.net
- [ ] Make the '>speces' to be the names of species used in 'tree.txt' file
- [ ] Make the '>resolution' to be the smallest size of your scaffolds (10kb)
- [ ] Make the '>numchr' to be the number of chromosomes in the closely-related reference species
*After adjusting, save as config.SFs*
```bash
# Directory that contains net files.
# The specified directory needs to have a sub-directory named as the reference species. For example, if the reference species is umd3 and the net files are placed in /aaa/bbb/ccc/umd3, then 'netdir' shoule be /aaa/bbb/ccc
>netdir
../Running_RACA/xiph.chain.net

# Directory that contains chain files.
# The specified directory needs to have a sub-directory named as the reference species. For example, if the reference species is umd3 and the chain files are placed in /aaa/bbb/ccc/umd3, then 'chaindir' shoule be /aaa/bbb/ccc
>chaindir
../Running_RACA/xiph.chain.net

# species-name tag (0: ref-species, 1: descendents, 2: outgroup)
>species
xiph 0
ppar 1
olat 2

# block resolution (bp)
>resolution
10000

# number of chromosomes in a ref-species (up to the X chromosome, excluding the Y chromosome)
>numchr
24
```

### 3b. Update and save a 'insertsize_sd.txt' file

- [ ] Updated the following tsv file with proper sample names and insert library sizes from the respective 'mapping_stats' files
- [ ] Save this file as /Running_RACA/insertsize_sd.txt

```bash
# Column1: insert library name
# Column2: insert library size from experiments (can be the same as Column3)
# Column3: mean size of insert library from read mapping (obtained from mapping_stats file)
# Column4: standard deviation (-) of insert library size (obtained from mapping_stats file)
# Column5: standard deviation (+) of insert library size (obtained from mapping_stats file)
# The positive and negative standard deviations could be the same
P26	2301	2301	132	132
P32	2300	2300	35	35
P34	2303	2303	114	114
P36	2303	2303	110	110
P37	2303	2303	111	111
P16	308	308	93	93
P25	367	367	99	99
P33	271	271	90	90
P35	337	337	100	100
P39	246	246	92	92
P26_15	15246	15246	154	154
P32_15	15299	15299	14	14
P34_15	15259	15259	143	143
P36_15	15262	15262	136	136
P37_15	15265	15265	133	133
```

### 3c. Update the following and save as 'params.txt' file
*This is where most of the settings for actually running RACA come from*
- [ ] Update 'INSERTLIBFILE' to point to your full path
- [ ] Update 'INSERTSIZETHR' to be a cut off in size between inserts of your paired end and mate paired reliable_adjs
- [ ] Update 'READMAPPINGDIR' to point to your full path
- [ ] Update 'READMAPPINGLIB' to point to your full path
- [ ] Change the 'NCPUS' to the number of threads you want RACA to use when it runs (I was running on 20 and it was taking 5-6 hrs)
- [ ] Update 'SCFSIZEFILE' to point to your full path
- [ ] Make sure your target names all start with "scaffold" or else change 'SCFPREFIX'
- [ ] Update 'SCFSEQFILE' to point to your full path
- [ ] Update 'OUTPUTDIR' to point to your full path
- [ ] Update 'TREEFILE' to point to your full path
- [ ] Update 'BENADJFILE' to point to your full path
- [ ] Update 'CONFIGSFSFILE' to point to your full path
- [ ] Update 'MAKESFSFILE' to point to your full path
- [ ] Save this as /Running_RACA/params.txt

```bash
# File that has the lengths of insert libraries and their means and standard deviations estimated from read mapping.
# Refer to the sample file 'insertsize_sd.txt'.
INSERTLIBFILE=../Running_RACA/insertsize_sd.txt

# Insert library size threshold for the normal directions of two end reads
# Size < INSERTSIZETHR : + -
# Otherwise : - +
# Due to the difference of library creation for short and long insert libraries
# If you think your insert libraries don't care about this, then use very large or small values to use the same criteria for read directions
INSERTSIZETHR=1000

# Input directory that has the paired-end read mapping data (mapping_file files produced from sam2RACA.pl)
# Refer to the file format by looking at the files in the TAreads directory.
# Current version only support that format.
# If you want to use any existing read alignment programs, you can simply convert the output from those alignment programs to the format that is supported by the current version of RACA.
READMAPPINGDIR=../Running_RACA/reads

# File that has the insert library name of each paired-end read mapping file in the $READMAPPING directory.
# Refer to the sample file 'data/readmapping_lib.txt'.
READMAPPINGLIB=../Running_RACA/readmapping_lib.txt

# The number of processes for parallel execution
NCPUS=20

# Size of target scaffolds
# Refer to the sample file 'panHod2.size'.
SCFSIZEFILE=../Running_RACA/Parae_female_10kb_assembly.size

# Prefix of target scaffold name
# If the name of a scaffold is "Scaffold1221", then SCFPREFIX should be Scaffold.
SCFPREFIX=scaffold

# Target scaffold sequences (this is the target_assembly.fa)
# This file contains all scaffold sequences
SCFSEQFILE=../Running_RACA/Parae_female_greater_10kb.fa

# Reference species
REFSPC=xiph

# Target species
TARSPC=ppar

# Window size for estimating paired-end read coverage threshold
WINDOWSIZE=1000

# Output directory
OUTPUTDIR=../RACA_Xiph_Parae_female/Running_RACA/Out_RACA

# Block resolution (bp)
RESOLUTION=10000

# The minimum percentage in a null distribution of P_ia(i,j) scores  
# that are obtained from entire scaffolds
# The actual P_ia(i,j) value in a null distribution that corresponds to
# MIN_INTRACOV_PERC is used as the cutoff threshold for P_ia(i,j)
# BEN FIND OUT WHAT THIS MEANS. LEAVING FOR NOW.
MIN_INTRACOV_PERC=5

# Sometimes, SF adjacencies only have comparative genomic information without
# paired-end read information because of the above MIN_INTRACOV_PERC threshold
# or long distance between two SFs.
# If this parameter is set (1), the SF adjacencies with both comparative
# genomic information and paired-end read information are used in the
# reconstruction
IGNORE_ADJS_WO_READS=0

# Newick tree file
# Refer to the sample file 'tree.txt'.
#
# Please append '@' symbol at the end of a target species name
#
TREEFILE=../Running_RACA/tree.txt

# Benchmark adjacency file
# Refer to the sample file 'reliable_adjs.txt'.
# If you don't have a benchmarking data for this file, then just specify an empty file. This will give an equal weight to the two components in the RACA's scoring function.
BENADJFILE=../Running_RACA/reliable_adjs.txt

# Config and make files for syntenic fragment construction
# Refer to the sample files 'config.SFs' and 'Makefile.SFs'.
# You need to change settings in the sample configuration file (config.SFs) according to your data
# You don't need to change anything in the Makefile.SFs.
CONFIGSFSFILE=../Running_RACA/config.SFs
MAKESFSFILE=../Running_RACA/Makefile.SFs
```

### 3d. Update the following and save as 'readmapping_lib.txt'
- [ ] This reflects the names of your mapping files (left column) and the names of your samples used in the 'insertsize_sd' file (right column)
- [ ] Save as '../Running_RACA/readmapping_lib.txt'
```bash
# Column1: file name
# Column2: insert library
P32_26_rc1_mapping_file	P26
P32_32_rc1_mapping_file	P32
P32_34_rc1_mapping_file	P34
P32_36_rc1_mapping_file	P36
P32_37_rc1_mapping_file	P37
P32_illum_16_mapping_file	P16
P32_illum_25_mapping_file	P25
P32_illum_33_mapping_file	P33
P32_illum_35_mapping_file	P35
P32_illum_39_mapping_file	P39
26_15kb_mapping_file	P26_15
32_15kb_mapping_file	P32_15
34_15kb_mapping_file	P34_15
36_15kb_mapping_file	P36_15
37_15kb_mapping_file	P37_15
```

### 3e. Make an empty 'reliable_adjs.txt' file
Save an empty file called '../Running_RACA/reliable_adjs.txt'
*This option is used by RACA if there are known adjustments to be made*

### 3f. Save a tree file to reflect node ages
Save as '../Running_RACA/tree.txt'
```bash
((xiph:42,ppar@:42):51,olat);
```

### 3g. Run RACA
```bash
perl Run_RACA.pl params.txt
```
