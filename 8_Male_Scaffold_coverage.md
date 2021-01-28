**Goals**
  To determine differences in alignment coverage by reads from each morph to the best *de novo* genome of each morph.

Approach:
  1) Align reads from each morph to the best *de novo* genome of each morph
  2) Determine scaffold coverage by each morph
  3) Compare scaffold coverage across morphs


# **1) Align reads from each morph to the best *de novo* genome of each morph**

## 1a) Index pseudohap files

```bash
for i in 32 34 36 37 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31
do
  bwa index P$i.1.fasta.gz
done
```

## 1b) Align reads to the best *de novo* genome of each morph

Parameters:
  bwa mem
    -t = number of threads
  samtools sort
    -n = Sort by read name
    -@4 = 4 threads
    -O BAM = Output in BAM format
    - = input and/or output is to the stdin
  samtools Fixmate
    -m = Add mate score tag
    -r Remove unmapped reads and secondary alignments
    -@4 = 4 threads
    -O BAM = Output in BAM format
    - = input and/or output is to the stdin

Best *de novo* genome of each morph:
  Immaculata - P09
  Parae - P04
  Melanzona - P01
**Note: from here on, repeat the entire pipeline for each of the best morph genomes- replace the word "morph" in script with relevant morph**

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 36 37 16 25 33 35 07 08 18 28 38 40
do
  screen -d -m bash -c "bwa mem -t 4 [best_morph_genome.fasta.gz] P$i\_for_paired.fastq.gz P$i\_rev_paired.fastq.gz | samtools sort -n -@4 -O BAM - | samtools fixmate -m -r -@4 -O BAM - [P$i\_to_morph.pseudohap.fixmate.bam]"
done
```

## 1c) Sort and mark duplicates

Parameters:
  samtools sort
    -@4 = 4 threads
    -O BAM = Output in BAM format
  samtools markdup
    -r Remove duplicate reads
    -S Mark supplementary alignments of duplicates as duplicates (slower). *Not sure why but everyone does this*
    -O BAM = Output in BAM format
    - = input and/or output is to the stdin

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 36 37 16 25 33 35 07 08 18 28 38 40
do
  screen -d -m bash -c "samtools sort -O bam -@4 [P$i\_to_morph.pseudohap.fixmate.bam] | samtools markdup  -@4 -r -S - [P$i\_to_morph.sorted.dedup.bam]"
done
```


# **2) Determine scaffold coverage by each morph**
## 2a) Run samtools depth with mapq60

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 37 16 25 35 07 18 28 38 40 36 33 08
do screen -d -m bash -c "
  samtools depth -aa -Q 60 [P$i\_to_morph.sorted.dedup.bam] | gzip > [P$i\_to_morph_q60.depth.txt.gz]"
done
```

## 2b) Get scaffold coverage by each sample

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 37 16 25 35 07 18 28 38 40 36 33 08
do
  zcat P$i\_to_P09_q60.depth.txt.gz | awk '{sum3[$1] += $3; count3[$1]++}; END{ for (id in sum3) { print id "\t" count3[id] "\t" sum3[id]/count3[id] } }' > [P$i\_to_morph_map60_raw.ScaffCov.txt]
done
```


## 2c) Correct raw average scaffold coverage by total scaffold coverage
*Note: This accounts for differences in sequencing library size across samples*

```bash
for i in 01 02 03 04 05 06 09 10 11 12 13 14 15 17 19 20 21 22 23 24 26 27 29 30 31 32 34 37 16 25 35 07 18 28 38 40 36 33 08
do
    avg=`awk '{sum += $3} END {print sum/NR}' [P$i\_to_morph_map60_raw.ScaffCov.txt]` ; awk -v avg=$avg '{print $0"\t"$3/avg}' [P$i\_to_morph_map60_raw.ScaffCov.txt] > [P$i\_to_morph_map60.ScaffCov.avg]
done
```


# **3) Compare scaffold coverage across morphs**
## 3a) Make directories for coverage by each morph
*Note: this is to use a script that will find coverage by all files in a directory*

```bash
mkdir coverage
cd coverage
mkdir blue red yellow mel males immac parae females nonimmac
cd ..
```

## 3b) Move individuals to their respective morphs
```bash
#Females
for i in 15 16 25 26 32 33 34 35 36 37 39
do
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/females/
done

#Melanzona
for i in 01 02 11 12 19 20 21 38 07 08 17 18 13 14 22 23 24
do
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/mel
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/nonimmac
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/nonparae
done

#Parae
for i in 03 04 27 28 29 30 31
do
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/parae
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/nonimmac
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/nonmel
done

#Immaculata
for i in 05 06 09 10 40
do
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/immac
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/nonparae
  cp [P$i\_to_morph_map60.ScaffCov.avg] coverage/nonmel
done
```

## 3c) Calculate scaffold coverage by each morph

```bash
python coverage_calc.py females/ [females_morph.q60_avg.txt]
python coverage_calc.py immac/ [immac_morph.q60_avg.txt]
python coverage_calc.py parae/ [parae_morph.q60_avg.txt]
python coverage_calc.py mel/ [mel_morph.q60_avg.txt]
python coverage_calc.py nonimmac/ [nonimmac_morph.q60_avg.txt]
python coverage_calc.py nonmel/ [nonmel_morph.q60_avg.txt]
python coverage_calc.py nonparae/ [nonparae_morph.q60_avg.txt]
```

## 3d) Determine fold change of scaffold coverage between morphs

```bash
python foldchange_coverage.py [females_morph.q60_avg.txt] [immac_morph.q60_avg.txt] [immac_vs_females_morph.q60_avg.txt]
python foldchange_coverage.py [females_morph.q60_avg.txt] [males_morph.q60_avg.txt] [males_vs_females_morph.q60_avg.txt]
python foldchange_coverage.py [nonimmac_morph.q60_avg.txt] [immac_morph.q60_avg.txt] [immac_vs_nonimmac_morph.q60_avg.txt]
```
