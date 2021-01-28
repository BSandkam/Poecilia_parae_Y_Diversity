**Goals**
  To identify scaffolds with morph unique sequence and prepare them for annotation.

Approach:
  1) Align morph-mers (morph specific Y-mers) to *de novo* genomes
  2) Find all scaffolds that contain at least 5 morph-mers


# **1) Align morph-mers (morph-specific Y-mers) to de novo genomes**
## 1a) For each morph generate file of morph-specific Ymers in just forward orientation
*Note: This uses the files of morph-mers generated in step 4 that start with "uniq_to"*
```bash
for FILE in uniq_to*
do
  awk '{print$1}' $FILE > $FILE.1col
done
```


## 1b) Make morph-mers into fasta format for alignment

```bash
for FILE in *.1col
do
  python3.5 seq_lines_to_fasta.py $FILE
done
```

## 1c) Align morph-mers to best *de novo* male genome of each morph

Parameters:
  bowtie2 (v2.4.1)
    -f = input files are fasta format
    -N = max # mismatches in seed alignment
    -L = length of seed substrings
    --end-to-end = entire read must align (no clipping)
    --all = report all alignments
    --threads = threads
    -x = index filename prefix
    -U = fasta of kmers to align

**Note: from here on, repeat the pipeline for each of the best morph genomes- replace the word "morph" in script with relevant morph and the word "sample" with the sample of that morph**
Best *de novo* genome of each morph:
  Immaculata - sample P09
  Parae - sample P04
  Melanzona - sample P01

```bash
bowtie2-build --threads 20 -f [sample.1.fasta.gz] [sample.1.fasta.bowtie2.index]
bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x [sample.1.fasta.bowtie2.index] -U [fasta_of_uniq_to_all_morph.1col.fa] | samtools sort -O bam -m 3G -o [sample.morphmers.bowtie.bam]
```


# **2) Find all scaffolds that contain at least 5 morph-mers**
## 2a) Convert BAM to BED
The BAM file then needs to be converted to BED for the next steps. The script `kmers-BAM2BED.py` can help with this. It should work with any BAM file but was designed specifically for this pipeline so it's likely not appropriate for other uses. I added an awk for $5==0 because without that I was still getting results from alignments with one mismatch.

Parameters:
  sort
    -S1G = buffer size for main memory
    -k1,1 -k2,2n = sort by column 1 first, then column2 (numerically)
  gzip
    -c = write to stdout

```bash
kmers-BAM2BED.py [sample.morphmers.bowtie.bam] | awk '($5==0){print$0}' | sort -S1G -k1,1 -k2,2n | gzip -c > [sample.morphmers.bowtie.bed.gz]
```

## 2b) Grab scaffolds with at least 5 morph-mers mapped

```bash
zcat [sample.morphmers.bowtie.bed.gz] | awk '{print$1}' | uniq -c | awk '($1>5){print$2}' | seqtk subseq [sample.1.fasta.gz] /dev/stdin \
| awk '{print $1}' \
| awk '{ if($1 ~ /^>/) {print $1"\_P01"} else {print $1} }' > [morph_scaffolds.fas]
```
