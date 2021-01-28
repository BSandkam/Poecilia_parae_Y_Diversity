**Goals**
  To identify morph specific sequence by comparing k-mer presence/absence across individuals, male morphs, sexes, and species

Approach:
  1) Identify all 31bp k-mers from all 29 individual de novo genomes
  2) Identify all 31bp k-mers from trimmed Illumina reads for all 40 individuals
  3) Remove female k-mers from each male, and male k-mers from each female
  4) Find Y-mers (male k-mers not present in females) that are present in every individual of each morph
  5) Find Y-mers unique to each morph (and present in every individual of that morph)
  6) Validate pipeline


*Note: working with the kmer files in this pipeline takes a lot of compute power (at times 600GB of RAM). This was all run on system with 48 cores and 1.5 Tb RAM and file sizes on disk at some stages was over 1Tb*

Samples (P*) by morph:
  *Red Melanzona*
    P01 P02 P11 P12 P19 P20 P21
  *Yellow Melanzona*
    P13 P14 P22 P23 P24
  *Blue Melanzona*
    P17
  *Parae morph*
    P03 P04 P27 P29 P30 P31
  *Immaculata*
    P05 P06 P09 P10
  *Female*
    P15 P26 P32 P34 P36 P37

# Step 1) Identify all 31bp k-mers from all 29 individual de novo genomes
*Note: Supernova can generate two output styles of de novo genome: megabubble (higher phasing but shorter scaffolds) or pseudohap (random flattening of meggabubbles to generate phased scaffolds that are longer but suffer from not all sequence being incorporated). To avoid kmer drop out during flattening we use the megabubble output for kmer comparisons.*


## 1a) Find all kmers in each indivdiual the Generate jellyfish files for each sample from megabubbles
Parameters
  jellyfish count (v2.2.3)
    -m = kmer length
    -C = find all canonical kmers (both forward and reverse strand of sequence) (used because we don't know if the same scaffolds are going the same direction across individuals)
    -s = initial hash size for running (how many kmers are stored in memory)
    -t = threads
    -o = output file
```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  jellyfish count -m 31 -C -s 8G -t 12 $i\_megbub.fasta.gz -o $i.meg.jf
done
```


## 1b) Generate jellyfish dump files (human readable) of the kmers
Parameters
  jellyfish dump (v2.2.3)
    -t = tab format so the number of occurances can be easily stripped to improve downstream computation
    -o = output dumped file
```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  jellyfish dump -t $i.meg.jf -o $i.meg.jf.dumped
done
```


## 1c) Remove the number of times each kmer occurs to improve downstream analyses

```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  grep -v ">" $i.meg.jf.dumped > $i.meg.just_kmers
done
```


## 1d) Use script Add_reverse_comp.py to add reverse complements to kmers

```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  python3.5 Add_reverse_comp.py -d $i.meg.just_kmers -o $i.meg.just_kmers_wcomp
done
```


# Step 2) Identify all 31bp k-mers from trimmed Illumina reads for all 40 individuals

## 2a) Find all kmers present in the trimmed Illumina reads of females and males
*Note: in addition to the 29 samples sequenced through the 10X linked-read pipeline, this also includes the 11 samples that were directly sequenced on the Illumina platform*

*Note: low count kmers from sequencing reads are often the result of sequencing errors. To avoid kmers from sequencing error we require kmers be present at least 4 times in an individual (sequencing was done to depth of ~40X so it is unliklely that true kmers were dropped)*

Parameters
  jellyfish count (v2.2.3)
    -m = kmer length
    -C = find all canonical kmers (both forward and reverse strand of sequence)
    -s = initial hash size for running (how many kmers are stored in memory)
    -t = threads
    -L 4 = Don't output k-mers with count <4
    -o = output file

```bash
for i in 01 03 05 09 11 13 15 19 21 23 26 29 31 34 37 02 04 06 10 12 14 17 20 22 24 27 30 32 36 16 25 33 35 39 07 08 18 28 38 40
do
  jellyfish count -m 31 -o P$i.fq.jelly -C -s 7G -t 12 -L 4 P$i\_for_paired.fastq P$i\_rev_paired.fastq
done
```


## 2b) Generate file of all k-mers present in any female and file of all k-mers present in any male

```bash
# Merging all females
jellyfish merge P15.fq.jelly P26.fq.jelly P32.fq.jelly P34.fq.jelly P36.fq.jelly P37.fq.jelly P16.fq.jelly P25.fq.jelly P33.fq.jelly P35.fq.jelly P39.fq.jelly -o females_fq_merged
#Generating female dump file to make human readable
jellyfish dump females_fq_merged -o females_fq_merged_dumped

# Merging all males
jellyfish merge P01.fq.jelly P03.fq.jelly P05.fq.jelly P09.fq.jelly P11.fq.jelly P13.fq.jelly P19.fq.jelly P21.fq.jelly P23.fq.jelly P29.fq.jelly P31.fq.jelly P02.fq.jelly P04.fq.jelly P06.fq.jelly P10.fq.jelly P12.fq.jelly P14.fq.jelly P17.fq.jelly P20.fq.jelly P22.fq.jelly P24.fq.jelly P27.fq.jelly P30.fq.jelly P07.fq.jelly P08.fq.jelly P18.fq.jelly P28.fq.jelly P38.fq.jelly P40.fq.jelly -o males_fq_merged

# Generating male dump file to make human readable
jellyfish dump males_fq_merged -o males_fq_merged_dumped
```


## 2c) Remove the number of times each kmer occurs to improve downstream analyses

```bash
sed '/>/d' females_fq_merged_dumped | sort -k 1,1 -S 40% --parallel=15 > females_fq_merged.noNum
sed '/>/d' males_fq_merged_dumped | sort -k 1,1 -S 40% --parallel=15 > males_fq_merged_dumped.noNum
```


## 2d) Use script Add_reverse_comp.py to add reverse complements to kmers

```bash
python3.5 Add_reverse_comp.py -d females_fq_merged.noNum -o females_fq_merged.noNum.revcomp.2col
python3.5 Add_reverse_comp.py -d males_fq_merged_dumped.noNum -o males_fq_merged.noNum.revcomp.2col
```



# Step 3) Remove female k-mers from each male, and male k-mers from each female
*Note: In theory, females should not have any k-mers that are not present in males. Therefore we run females through the same pipeline as males to provide false positive rate of male specific k-mer (Y-mer) identification*

## 3a) Remove k-mers directly generated from Illumina reads
```bash
# Remove female Illumina k-mers from male genome k-mers
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P17 P19 P20 P21 P22 P23 P24 P27 P29 P30 P31 P15 P26 P32 P34 P36 P37
do
  awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' females_fq_merged.noNum.revcomp.2col $i.meg.just_kmers_wcomp > $i.meg.kmers.not_in_female_fq
done

# Remove male Illumina k-mers from female genome k-mers
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P17 P19 P20 P21 P22 P23 P24 P27 P29 P30 P31 P15 P26 P32 P34 P36 P37
do
  awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.noNum.revcomp.2col $i.meg.just_kmers_wcomp > $i.meg.kmers.not_in_male_fq
done
```


## 3b) Remove k-mers generated from de novo genomes
*Note: this step is done to remove any true k-mers that may have not been present >3 times in the fastq files due to stocastic low coverage of a particular region*
```bash
# Join female megabubble kmers not in female fastq kmers
cat P32.meg.kmers.not_in_female_fq P34.meg.kmers.not_in_female_fq P36.meg.kmers.not_in_female_fq P37.meg.kmers.not_in_female_fq P15.meg.kmers.not_in_female_fq P26.meg.kmers.not_in_female_fq > all_female_not_in_female_fq

# Sort female megabubble kmers not in female fastq kmers
sort all_female_not_in_female_fq | uniq > all_female_not_in_female_fq_uniq

# Remove female megabubble kmers not in female fastq kmers from male megabubble kmers
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P17 P19 P20 P21 P22 P23 P24 P27 P29 P30 P31
do
  awk 'NR==FNR{a[$1];a[$2];next}!($1 in a){print $0}' all_female_not_in_female_fq_uniq $i.meg.kmers.not_in_female_fq > $i.meg.kmers.not_in_female
done


## Do the inverse of males and females to get false positive rate

# Join male megabubble kmers not in male fastq kmers
cat P01.meg.kmers.not_in_male_fq P02.meg.kmers.not_in_male_fq P03.meg.kmers.not_in_male_fq P04.meg.kmers.not_in_male_fq P05.meg.kmers.not_in_male_fq P06.meg.kmers.not_in_male_fq P09.meg.kmers.not_in_male_fq P10.meg.kmers.not_in_male_fq P11.meg.kmers.not_in_male_fq P12.meg.kmers.not_in_male_fq P13.meg.kmers.not_in_male_fq P14.meg.kmers.not_in_male_fq P17.meg.kmers.not_in_male_fq P19.meg.kmers.not_in_male_fq P20.meg.kmers.not_in_male_fq P21.meg.kmers.not_in_male_fq P22.meg.kmers.not_in_male_fq P23.meg.kmers.not_in_male_fq P24.meg.kmers.not_in_male_fq P27.meg.kmers.not_in_male_fq P29.meg.kmers.not_in_male_fq P30.meg.kmers.not_in_male_fq P31.meg.kmers.not_in_male_fq > all_male_not_in_male_fq

# Sort male megabubble kmers not in male fastq kmers
sort all_male_not_in_male_fq | uniq > all_male_not_in_male_fq_uniq

# Remove male megabubble kmers not in male fastq kmers from female megabubble kmers
for i in P15 P26 P32 P34 P36 P37
do
  awk 'NR==FNR{a[$1];a[$2];next}!($1 in a){print $0}' all_male_not_in_female_fq_uniq $i.meg.kmers.not_in_male_fq > $i.meg.kmers.not_in_male
done
```


# Step 4) Find Y-mers (male k-mers not present in females) that are present in every individual of each morph

## 4a) Find Y-mers present in every male
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P02.meg.kmers.not_in_female P02.meg.kmers.not_in_female > in_01_02
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P03.meg.kmers.not_in_female in_01_02 > in_01_02_03
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P04.meg.kmers.not_in_female in_01_02_03 > in_01_02_03_04
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P05.meg.kmers.not_in_female in_01_02_03_04 > in_01_02_03_04_05
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P06.meg.kmers.not_in_female in_01_02_03_04_05 > in_01_02_03_04_05_06
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P09.meg.kmers.not_in_female in_01_02_03_04_05_06 > in_01_02_03_04_05_06_09
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P10.meg.kmers.not_in_female in_01_02_03_04_05_06_09 > in_01_02_03_04_05_06_09_10
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P11.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10 > in_01_02_03_04_05_06_09_10_11
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P12.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11 > in_01_02_03_04_05_06_09_10_11_12
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P13.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12 > in_01_02_03_04_05_06_09_10_11_12_13
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P14.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13 > in_01_02_03_04_05_06_09_10_11_12_13_14
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P17.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14 > in_01_02_03_04_05_06_09_10_11_12_13_14_17
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P19.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P20.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P21.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P22.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P23.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P24.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P27.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24_27
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P29.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24_27 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24_27_29
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P30.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24_27_29 > in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24_27_29_30
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P31.meg.kmers.not_in_female in_01_02_03_04_05_06_09_10_11_12_13_14_17_19_20_21_22_23_24_27_29_30 > in_all_males
```


## 4b) Find Y-mers present in every Red melanzona
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P01.meg.kmers.not_in_female P02.meg.kmers.not_in_female > red_01_02
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P11.meg.kmers.not_in_female red_01_02 > red_01_02_11
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P12.meg.kmers.not_in_female red_01_02_11 > red_01_02_11_12
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P20.meg.kmers.not_in_female red_01_02_11_12 > red_01_02_11_12_20
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P21.meg.kmers.not_in_female red_01_02_11_12_20 > red_01_02_11_12_20_21
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P19.meg.kmers.not_in_female red_01_02_11_12_20_21 > in_all_red
```


## 4c) Find Y-mers present in every Yellow melanzona
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P13.meg.kmers.not_in_female P14.meg.kmers.not_in_female > yel_13_14
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P22.meg.kmers.not_in_female yel_13_14 > yel_13_14_22
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P23.meg.kmers.not_in_female yel_13_14_22 > yel_13_14_22_23
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P24.meg.kmers.not_in_female yel_13_14_22_23 > in_all_yellow
```


## 4d) Find Y-mers present in every Parae morph
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P03.meg.kmers.not_in_female P04.meg.kmers.not_in_female > parae_03_04
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P27.meg.kmers.not_in_female parae_03_04 > parae_03_04_27
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P29.meg.kmers.not_in_female parae_03_04_27 > parae_03_04_27_29
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P30.meg.kmers.not_in_female parae_03_04_27_29 > parae_03_04_27_29_30
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P31.meg.kmers.not_in_female parae_03_04_27_29_30 > in_all_parae_morph
```


## 4e) Find Y-mers present in every Immaculata morph
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P05.meg.kmers.not_in_female P06.meg.kmers.not_in_female > immac_05_06
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P09.meg.kmers.not_in_female immac_05_06 > immac_05_06_09
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P10.meg.kmers.not_in_female immac_05_06_09 > in_all_immaculata
```


## 4f) Find Y-mers present in every melanzona (red, blue, and yellow)
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' in_all_red in_all_yellow > in_all_RED_and_YELLOW
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' in_all_RED_and_YELLOW P17.meg.kmers.not_in_female > in_all_melanzona
```


## 4g) Find Y-mers present in every female
*Note: this is done as confirmation of the pipeline since any female unique k-mers should be individual SNPs and therefore not shared between all females. As expected we found 0 female unique kmers that were present in every female.*
```bash
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P15.meg.kmers.not_in_male P26.meg.kmers.not_in_male > inF_15_26
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P32.meg.kmers.not_in_male inF_15_26 > inF_15_26_32
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P34.meg.kmers.not_in_male inF_15_26_32 > inF_15_26_32_34
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P36.meg.kmers.not_in_male inF_15_26_32_34 > inF_15_26_32_34_36
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P37.meg.kmers.not_in_male inF_15_26_32_34_36 > in_all_females
```



# Step 5) Find Y-mers unique to each morph - and present in every individual of that morph
*This step finds what Y-mers that are present in every individual of a morph are NOT present an ANY of the other individuals of the other morphs*
*Note: these all work by combining Y-mers from all individuals that are not of the designated morph into a file called 'anti' then take all Ymers in the anti-file out of the file of Y-mers that are present in every individual of that morph.*

## 5a) Find Red unique Y-mers
```bash
for i in P03 P04 P05 P06 P09 P10 P13 P14 P15 P17 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_red
done
sort -u anti_red > anti_red_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_red_uniq in_all_red > uniq_to_red
```


## 5b) Find Yellow unique Y-mers
```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P15 P17 P19 P20 P21 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_yellow
done
sort -u anti_yellow > anti_yellow_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_yellow_uniq in_all_yellow > uniq_to_yellow
```


## 5c) Find Blue unique Y-mers
```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P15 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_blue
done
sort -u anti_blue > anti_blue_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_blue_uniq P17.meg.kmers.not_in_female > uniq_to_blue
```


## 5d) Find Parae morph unique Y-mers
```bash
for i in P01 P02 P05 P06 P09 P10 P11 P12 P13 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_parae
done
sort -u anti_parae > anti_parae_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_parae_uniq in_all_parae_morph > uniq_to_parae_morph
```


## 5e) Find Immaculata unique Y-mers
```bash
for i in P01 P02 P03 P04 P11 P12 P13 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_immac
done
sort -u anti_immac > anti_immac_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_immac_uniq in_all_immaculata > uniq_to_immaculata
```


## 5f) Find Melanzona unique Y-mers
```bash
for i in P03 P04 P05 P06 P09 P10 P15 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_mel
done
sort -u anti_mel > anti_mel_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_mel_uniq in_all_melanzona > uniq_to_melanzona
```



# Step 6) Validate pipeline

## 6a) Validate morph unique Y-mers are not due to chance based on sample size
Validations to run:
- 4 individuals (no blue) : P01 P13 P03 P05
- 5 individuals (one of each morph) : P01 P13 P17 P03 P05
- 9 individuals (two of each morph but the one blue) : P02 P11 P14 P22 P17 P04 P27 P09 P10
- 13 individuals (three of each morph but the one blue) : P12 P19 P20 P22 P23 P24 P17 P29 P30 P31 P06 P09 P10

### 6a.1) Find Ymers unique to 4 random individuals: P01 P13 P03 P05
```bash
# Find Ymers in all 4 individuals
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P01.meg.kmers.not_in_female P13.meg.kmers.not_in_female > valid_01_13
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P03.meg.kmers.not_in_female valid_01_13 > valid_01_13_03
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P05.meg.kmers.not_in_female valid_01_13_03 > valid_4indiv

# Find the uniqueness of N=4 validation
for i in P02 P04 P06 P09 P10 P11 P12 P14 P15 P17 P19 P20 P21 P22 P23 P24 P26 P27 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_valid4
done
sort -u anti_valid4 > anti_valid4_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_valid4_uniq valid_4indiv > uniq_valid_4indiv
```


### 6a.2) Find Ymers unique to 5 random individuals: P12 P23 P17 P29 P06
```bash
# Finding Ymers in all 5 individuals
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P12.meg.kmers.not_in_female P23.meg.kmers.not_in_female > valid_12_23
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P17.meg.kmers.not_in_female valid_12_23 > valid_12_23_17
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P29.meg.kmers.not_in_female valid_12_23_17 > valid_12_23_17_29
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P06.meg.kmers.not_in_female valid_12_23_17_29 > valid_5indiv

# Finding the uniqueness of N=5 validation
for i in P01 P13 P03 P05 P02 P04 P09 P10 P11 P14 P15 P19 P20 P21 P22 P24 P26 P27 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_valid5
done
sort -u anti_valid5 > anti_valid5_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_valid5_uniq valid_5indiv > uniq_valid_5indiv
```


### 6a.3) Find Ymers unique to 9 random individuals: P02 P11 P14 P22 P17 P04 P27 P09 P10
```bash
# Finding Ymers in all 9 individuals
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P02.meg.kmers.not_in_female P11.meg.kmers.not_in_female > valid_02_11
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P14.meg.kmers.not_in_female valid_02_11 > valid_02_11_14
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P22.meg.kmers.not_in_female valid_02_11_14 > valid_02_11_14_22
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P17.meg.kmers.not_in_female valid_02_11_14_22 > valid_02_11_14_22_17
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P04.meg.kmers.not_in_female valid_02_11_14_22_17 > valid_02_11_14_22_17_04
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P27.meg.kmers.not_in_female valid_02_11_14_22_17_04 > valid_02_11_14_22_17_04_27
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P09.meg.kmers.not_in_female valid_02_11_14_22_17_04_27 > valid_02_11_14_22_17_04_27_09
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P10.meg.kmers.not_in_female valid_02_11_14_22_17_04_27_09 > valid_9indiv

# Finding the uniqueness of N=9 validation
for i in P01 P03 P05 P06 P12 P13 P15 P19 P20 P21 P23 P24 P26 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_valid9
done
sort -u anti_valid9 > anti_valid9_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_valid9_uniq valid_9indiv > uniq_valid_9indiv
```


### 6a.4) Find Ymers unique to 13 random individuals: P12 P19 P20 P22 P23 P24 P17 P29 P30 P31 P06 P09 P10
```bash
# Finding Ymers in all 13 individuals
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P12.meg.kmers.not_in_female P19.meg.kmers.not_in_female > valid_12_19
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P20.meg.kmers.not_in_female valid_12_19 > valid_12_19_20
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P22.meg.kmers.not_in_female valid_12_19_20 > valid_12_19_20_22
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P23.meg.kmers.not_in_female valid_12_19_20_22 > valid_12_19_20_22_23
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P24.meg.kmers.not_in_female valid_12_19_20_22_23 > valid_12_19_20_22_23_24
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P17.meg.kmers.not_in_female valid_12_19_20_22_23_24 > valid_12_19_20_22_23_24_17
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P29.meg.kmers.not_in_female valid_12_19_20_22_23_24_17 > valid_12_19_20_22_23_24_17_29
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P30.meg.kmers.not_in_female valid_12_19_20_22_23_24_17_29 > valid_12_19_20_22_23_24_17_29_30
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P31.meg.kmers.not_in_female valid_12_19_20_22_23_24_17_29_30 > valid_12_19_20_22_23_24_17_29_30_31
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P06.meg.kmers.not_in_female valid_12_19_20_22_23_24_17_29_30_31 > valid_12_19_20_22_23_24_17_29_30_31_06
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P09.meg.kmers.not_in_female valid_12_19_20_22_23_24_17_29_30_31_06 > valid_12_19_20_22_23_24_17_29_30_31_06_09
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' P10.meg.kmers.not_in_female valid_12_19_20_22_23_24_17_29_30_31_06_09 > valid_13indiv

# Finding the uniqueness of N=13 validation
for i in P01 P02 P03 P04 P05 P11 P13 P14 P15 P21 P26 P27 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_valid13
done
sort -u anti_valid13 > anti_valid13_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_valid13_uniq valid_13indiv > uniq_valid_13indiv
```

### 6a.5) Find Ymers unique to overlap of random subsamples from above (N=4, N=5, N=9)
```bash
# Find Ymers present in all individuals of subsets
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' valid_4indiv valid_5indiv > valid_in_all_4_and_5
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' valid_4indiv valid_9indiv > valid_in_all_4_and_9
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' valid_5indiv valid_9indiv > valid_in_all_5_and_9
awk 'NR==FNR{a[$1];next}($1 in a) || ($2 in a){print $0}' valid_in_all_4_and_5 valid_9indiv > valid_in_all_4_5_and_9

# Find uniqueness of Ymers in random 4 and 5 individual groups
for i in P02 P04 P09 P10 P11 P14 P15 P19 P20 P21 P22 P24 P26 P27 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_4_and_5
done
sort -u anti_4_and_5 > anti_4_and_5_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_4_and_5_uniq valid_in_all_4_and_5 > uniq_to_valid_in_all_4_and_5

# Find uniqueness of Ymers in random 4 and 9 individual groups
for i in P06 P12 P15 P19 P20 P21 P23 P24 P26 P29 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_4_and_9
done
sort -u anti_4_and_9 > anti_4_and_9_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_4_and_9_uniq valid_in_all_4_and_9 > uniq_to_valid_in_all_4_and_9

# Find uniqueness of Ymers in random 5 and 9 individual groups
for i in P01 P03 P05 P13 P15 P19 P20 P21 P24 P26 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_5_and_9
done
sort -u anti_5_and_9 > anti_5_and_9_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_5_and_9_uniq valid_in_all_5_and_9 > uniq_to_valid_in_all_5_and_9

# Find uniqueness of Ymers in random 4, 5 and 9 individual groups
for i in P15 P19 P20 P21 P24 P26 P30 P31 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_4_5_and_9
done
sort -u anti_4_5_and_9 > anti_4_5_and_9_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_4_5_and_9_uniq valid_in_all_4_5_and_9 > uniq_to_valid_in_all_4_5_and_9
```


## 6.b) Validate to make sure that the Y-mers that are present in all males are not present in females
*Note: since we previously removed all female kmers we shouldn't see any differences between the in_all_males and the uniq_to_all_males and indeed we found no difference*
**Y-mers in all males are not present in female k-mers**
```bash
# Remove all kmers from the female samples from the Y-mers present in every male
for i in P15 P26 P32 P34 P36 P37
do
  cat $i.meg.kmers.not_in_female >> anti_males
done
sort -u anti_males > anti_males_uniq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' anti_males_uniq in_all_males > uniq_to_all_males
```


## 6.c) Make sure the morph unique Y-mers (morph-mers) don't occur in the de novo female genomes
*Approach: attempt to align the morph-mers to the female genomes with no mismatches or spaces*
**Morph-mers did not align to female genomes**

Parameters:
  bowtie2
    -f = input files are fasta format
    -N = max # mismatches in seed alignment
    -L = length of seed substrings
    --end-to-end = entire read must align (no clipping)
    --all = report all alignments
    --threads = threads
    -x = index filename prefix
    -U = fasta of kmers to align

```bash
# Make morph-mers into fasta to be aligned
for i in uniq_to_*
do
  awk '{print$1}' $i > $i.1col
done
for i in *.1col
do
  python3 seq_lines_to_fasta.py $i
done

# Map morph-mers to female genomes
for i in 15 26 32 34 36 37 01
do
  bowtie2-build --threads 20 -f P$i.1.fasta.gz P$i.1.fasta.bowtie2.index
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_all_males.1col.fa | samtools sort -O bam -m 3G -o P$i.ymers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_immaculata.1col.fa | samtools sort -O bam -m 3G -o P$i.immac-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_parae_morph.1col.fa | samtools sort -O bam -m 3G -o P$i.parae-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_melanzona.1col.fa | samtools sort -O bam -m 3G -o P$i.mel-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_red.1col.fa | samtools sort -O bam -m 3G -o P$i.red-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_blue.1col.fa | samtools sort -O bam -m 3G -o P$i.blue-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_yellow.1col.fa | samtools sort -O bam -m 3G -o P$i.yellow-mers.bowtie.bam
done

# Find if any morph-mers mapped perfectly
for i in 15 26 32 34 36 37 01
do
  kmers-BAM2BED.py P$i.ymers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.ymers
  kmers-BAM2BED.py P$i.immac-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.immac-mers
  kmers-BAM2BED.py P$i.parae-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.parae-mers
  kmers-BAM2BED.py P$i.mel-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.mel-mers
  kmers-BAM2BED.py P$i.red-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.red-mers
  kmers-BAM2BED.py P$i.yellow-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.yellow-mers
  kmers-BAM2BED.py P$i.blue-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.blue-mers
done
```


## 6.d) Make sure the morph-mers don't occur in genomes of males of other morphs
*Approach: align morph mers to the males of the other morphs. There should be none aligned*
Parameters:
  bowtie2
    -f = input files are fasta format
    -N = max # mismatches in seed alignment
    -L = length of seed substrings
    --end-to-end = entire read must align (no clipping)
    --all = report all alignments
    --threads = threads
    -x = index filename prefix
    -U = fasta of kmers to align (created in previous step)

```bash
# Align morph mers
for i in 01 02 03 04 05 06 09 10 11 12 13 14 17 19 20 21 22 23 24 27 29 30 31
do
  bowtie2-build --threads 20 -f P$i.1.fasta.gz P$i.1.fasta.bowtie2.index
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_all_males.1col.fa | samtools sort -O bam -m 3G -o P$i.ymers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_immaculata.1col.fa | samtools sort -O bam -m 3G -o P$i.immac-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_parae_morph.1col.fa | samtools sort -O bam -m 3G -o P$i.parae-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_melanzona.1col.fa | samtools sort -O bam -m 3G -o P$i.mel-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_red.1col.fa | samtools sort -O bam -m 3G -o P$i.red-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_blue.1col.fa | samtools sort -O bam -m 3G -o P$i.blue-mers.bowtie.bam
  bowtie2 -f -N 0 -L 31 --end-to-end --all --threads 20 -x P$i.1.fasta.bowtie2.index -U fasta_of_uniq_to_yellow.1col.fa | samtools sort -O bam -m 3G -o P$i.yellow-mers.bowtie.bam
done

# Find how many morph-mers aligned
for i in 01 02 03 04 05 06 09 10 11 12 13 14 17 19 20 21 22 23 24 27 29 30 31
do
  kmers-BAM2BED.py P$i.ymers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.ymers
  kmers-BAM2BED.py P$i.immac-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.immac-mers
  kmers-BAM2BED.py P$i.parae-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.parae-mers
  kmers-BAM2BED.py P$i.mel-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.mel-mers
  kmers-BAM2BED.py P$i.red-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.red-mers
  kmers-BAM2BED.py P$i.yellow-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.yellow-mers
  kmers-BAM2BED.py P$i.blue-mers.bowtie.bam | awk '($5==0){print$0}' | wc -l > P$i.num.blue-mers
done

#Check to be sure none aligned - all these files should all be empty
for i in 03 04 05 06 09 10 13 14 17 22 23 24 27 29 30 31
do
  cat P$i.num.red-mers
done

for i in 01 02 03 04 05 06 09 10 11 12 17 19 20 21 27 29 30 31
do
  cat P$i.num.yellow-mers
done

for i in 01 02 03 04 05 06 09 10 11 12 13 14 19 20 21 22 23 24 27 29 30 31
do
  cat P$i.num.blue-mers
done

for i in 01 02 05 06 09 10 11 12 13 14 17 19 20 21 22 23 24
do
  cat P$i.num.parae-mers
done

for i in 01 02 03 04 11 12 13 14 17 19 20 21 22 23 24 27 29 30 31
do
  cat P$i.num.immac-mers
done

for i in 03 04 05 06 09 10 27 29 30 31
do
  cat P$i.num.mel-mers
done
```
