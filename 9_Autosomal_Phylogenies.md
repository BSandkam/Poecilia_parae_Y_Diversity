**Goals**
  To generate a phylogeny built on known autosomal sequence to ensure the pattern of Y diversity is not cryptic subpopulations

Approach
  1) Align reads from all samples to best de novo female genome (sample P32) to generate BAM files
  2) Find the longest female scaffold from each autosome (identified by RACA), and double check male:female was 1:1 in previous analysis (to further ensure it is autosomal sequence)
  3) Generate consensus sequences for each individual (N=40) from the longest scaffold on each chromosome
  4) Create alignments of consensus from all individuals for the longest scaffold on each chromosome
  5) Build phylogenies for each chromosome


# **Step 1) Align reads from all samples to best de novo female genome (sample P32) to generate BAM files**

## 1a) Align reads to P32 pseudohap

Parameters:
  bwa mem
    -t = number of threads
  samtools sort
    -n = Sort by read name
    -@1 = 1 threads
    -O BAM = Output in BAM format
    - = input and/or output is to the stdin
  samtools Fixmate
    -m = Add mate score tag
    -r Remove unmapped reads and secondary alignments
    -@1 = 1 threads
    -O BAM = Output in BAM format
    - = input and/or output is to the stdin

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 36 37 25 33 35 39 07 08 18 40 16 38 28
do
  screen -d -m bash -c "bwa mem -t 4 P32_greater_10kb.fa P$i\_for_paired.fastq.gz P$i\_rev_paired.fastq.gz | samtools sort -n -@4 -O BAM - | samtools fixmate -m -r -@4 -O BAM - P$i\_to_P32.fixmate.bam"
done
```


## 1b) Sort and mark duplicates

Parameters:
  samtools sort
    -@1 = 1 thread
    -O BAM = Output in BAM format
  samtools markdup
    -r Remove duplicate reads
    -S Mark supplementary alignments of duplicates as duplicates (slower).
    -O BAM = Output in BAM format
    - = input and/or output is to the stdin

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 36 37 25 33 35 39 07 08 18 40 16 38 28
do
  screen -d -m bash -c "samtools sort -O bam -@4 P$i\_to_P32.fixmate.bam | samtools markdup  -@4 -r -S - P$i\_to_P32.sorted.dedup.bam"
done
```


# **Step 2) Find the longest female scaffold from each autosome (identified by RACA), and double check male:female was 1:1 in previous analysis (to further ensure it is autosomal sequence)**

## 2a) Find the longest female scaffold from each autosome (identified by RACA)
*Note: this uses the segments file generated in the output of the RACA pipeline*
*Note: change the awk'($1~1.) so the ~1. is the chromosome of focus*

```bash
grep "RACA." rec_chrs.ppar.segments.refined.txt | awk '{print $1}' | awk -F '[.:-]' '{print $2  " " $4-$3 " " $3 "-" $4}'| sort -k 2 -nr | awk '($1~1.){print}' | head
```


## 2b) Find average coverage of each of these scaffolds to verify coverage is normal for both males and females and male:female is ~ 1:1 (further confirming scaffold is from autosome)
*Note: this uses the average-coverage-fold-change file generated during the female coverage pipeline*
*Note: replace scaffold number for each scaffold identified above*
```bash
grep "scaffold163148" male_female_avg_coverage_fold_change.txt | awk -F"," '{print $1 "\t" "MFlog=" "\t" $2 "\t" "Male=" "\t" $3 "\t" "Female=" "\t" $5}'
```


# **Step 3) Generate consensus sequences for each individual (N=40) from the longest scaffold on each chromosome**

## 3a) Use BCFtools to generate consensus sequences for each scaffold
Parameters
  bcftools (v1.11-35-g8a744dd)
    mpileup
      -Ou = output uncompressed BCF
      -f = indexed reference file
    call
      -mv = default calling method with output only variant sites
      -Ob = output compressed BCF format
      -o = output file
```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 36 37 25 33 35 39 07 08 18 40 16 38 28
do
  screen -d -m bash -c "
  # call variants
  bcftools mpileup -Ou -f P32_greater_10kb.fa P$i\_to_P32.sorted.dedup.bam | bcftools call -mv -Oz -o P$i\_to_P32_calls.vcf.gz

  bcftools index P$i\_to_P32_calls.vcf.gz

  samtools faidx -r scaffolds2.bed P32_greater_10kb.fa | bcftools consensus -I P$i\_to_P32_calls.vcf.gz -o P$i\_to_P32_consensus_IUPAC.fa"
done
```

## 3b) Generate fasta files for each chromosome with sample name as header instead of scaffold number

```bash
for i in {1..24}
do
  mkdir chr$i
done

#Index all consensus sequences
for i in P01 P02 P03 P04 P05 P06 P07 P09 P10 P11 P12 P13 P14 P15 P17 P18 P19 P20 P21 P22 P23 P24 P25 P27 P26 P30 P31 P32 P34 P29 P08 P16 P28 P33 P35 P36 P37 P38 P39 P40
do
  samtools faidx $i\_to_P32_consensus_IUPAC.fa
done

#Note: Run the next step for each chromosome, adjust directory and scaffold - Chromosome 1 as example

# Chr1
for i in P01 P02 P03 P04 P05 P06 P07 P09 P10 P11 P12 P13 P14 P15 P17 P18 P19 P20 P21 P22 P23 P24 P25 P27 P26 P30 P31 P32 P34 P29 P08 P16 P28 P33 P35 P36 P37 P38 P39 P40
do
  samtools faidx $i\_to_P32_consensus_IUPAC.fa scaffold163148 | sed "1s/.*/>${i}/"  > chr1/$i.fa
done
```


# **Step 4) Create alignments of consensus from all individuals for the longest scaffold on each chromosome**
Each directory was imported to Geneious Prime as a sequence list.
Alignments of each scaffold were then made with the MAFFT plug-in for Geneious.


# **Step 5) Build phylogenies for each chromosome**
Phylogenies of each scaffold were made with the FastTree plug-in for Geneious Prime.
