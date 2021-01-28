**Goals:**
  Determine the false positive rate of Ymer identification. Run the female samples through the same narrowing pipeline used for males.

Overview of Steps:
  (1) Find all female kmers from the megabubble outputs (done as part of previous script)
  (2) Find all kmers from the fastq files of all males
  (3) Remove all male fastq kmers from the female megabubble kmers
  (4) Find all male megabubble kmers
  (5) Remove all male megabubble kmers from female megabubble kmers
  (6) Generate a matrix based on presence/absence of each remaining female megabubble kmer in each individual
  (7) Find all femmers that are present in different numbers of females

# (2) Find all kmers from the fastq files of all males

## (2a) Making Male fastq kmers
 10x Reads stripped and trimmed
```bash
for i in 01 03 05 09 11 13 19 21 23 29 31 02 04 06 10 12 14 17 20 22 24 27 30 07 08 18 28 38 40
do
  jellyfish count -m 31 -o P$i.fq.jelly -C -s 7G -t 12 -L 4 <(zcat P$i\_for_paired.fastq.gz) <(zcat P$i\_rev_paired.fastq.gz)
done
```

## (2b) Merge and dump all MALE jelly files

```bash
jellyfish merge P01.fq.jelly P03.fq.jelly P05.fq.jelly P09.fq.jelly P11.fq.jelly P13.fq.jelly P19.fq.jelly P21.fq.jelly P23.fq.jelly P29.fq.jelly P31.fq.jelly P02.fq.jelly P04.fq.jelly P06.fq.jelly P10.fq.jelly P12.fq.jelly P14.fq.jelly P17.fq.jelly P20.fq.jelly P22.fq.jelly P24.fq.jelly P27.fq.jelly P30.fq.jelly P07.fq.jelly P08.fq.jelly P18.fq.jelly P28.fq.jelly P38.fq.jelly P40.fq.jelly -o males_fq_merged

jellyfish dump males_fq_merged -o males_fq_merged_dumped
```

## (2c) Make dump file just kmers

```bash
grep -v ">" males_fq_merged_dumped | sort -k 1,1 -S 40% --parallel=15 > males_fq_merged.just_kmers
```

## (2d) Add reverse complements MALE kmers

```bash
python3.5 Add_reverse_comp.py -d males_fq_merged.just_kmers -o males_fq_merged.just_kmers_wcomp
```


# (3) Remove all male fastq kmers from the female megabubble kmers

```bash
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.just_kmers_wcomp <(zcat /home/sandkam/kmers_April2020/P32.meg.just_kmers_wcomp.gz) > female32.meg.kmers_not_in_male_fq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.just_kmers_wcomp <(zcat /home/sandkam/kmers_April2020/P34.meg.just_kmers_wcomp.gz) > female34.meg.kmers_not_in_male_fq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.just_kmers_wcomp <(zcat /home/sandkam/kmers_April2020/P36.meg.just_kmers_wcomp.gz) > female36.meg.kmers_not_in_male_fq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.just_kmers_wcomp <(zcat /home/sandkam/kmers_April2020/P37.meg.just_kmers_wcomp.gz) > female37.meg.kmers_not_in_male_fq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.just_kmers_wcomp <(zcat /home/sandkam/kmers_April2020/P26.meg.just_kmers_wcomp.gz) > female26.meg.kmers_not_in_male_fq
awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' males_fq_merged.just_kmers_wcomp <(zcat /home/sandkam/kmers_April2020/P15.meg.just_kmers_wcomp.gz) > female15.meg.kmers_not_in_male_fq
cat female32.meg.kmers_not_in_male_fq female34.meg.kmers_not_in_male_fq female36.meg.kmers_not_in_male_fq female37.meg.kmers_not_in_male_fq female26.meg.kmers_not_in_male_fq female15.meg.kmers_not_in_male_fq | uniq > all_uniq_female.meg.kmers_not_in_male_fq
```


# (4) Find all male megabubble kmers

```bash
for i in P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P17 P19 P20 P21 P22 P23 P24 P27 P29 P30 P31
do
  zcat $i.meg.just_kmers_wcomp.gz >> all_male_meg_kmers
done
sort all_male_meg_kmers | uniq > all_UNIQ_male_meg_kmers
```


# (5) Remove all male megabubble kmers from female megabubble kmers

```bash
cd /home/sandkam/kmer_from_fastq/males
for i in 15 26 34 37 32 36
do
  cat female$i.meg.kmers_not_in_male_fq >> all_femmers
done
sort all_femmers | uniq > all_uniq_femmers

awk 'NR==FNR{a[$1];next}!($1 in a) && !($2 in a){print $0}' /home/sandkam/kmers_April2020/all_UNIQ_male_meg_kmers all_uniq_femmers > all_uniq_female.meg.kmers_not_in_male_fq_not_in_male_meg
awk '{print $1}' all_uniq_female.meg.kmers_not_in_male_fq_not_in_male_meg > all_uniq_femmers_1col
awk '{print $2}' all_uniq_female.meg.kmers_not_in_male_fq_not_in_male_meg >> all_uniq_femmers_1col
sort all_uniq_femmers_1col | uniq > all_uniq_femmers_1col_uniq
```


# (6) Generate a matrix based on presence/absence of each remaining female megabubble kmer in each individual

## (6a) Generate file of presence absence of femmer for each individual

```bash
cd /home/sandkam/kmer_from_fastq/males
for i in 15 26 34 37 32 36
do
awk 'NR==FNR{a[$1];a[$2];next} {
  if (($1 in a) || ($2 in a))
  print "1";
else
  print "0";
}' female$i.meg.kmers_not_in_male_fq all_uniq_femmers_1col_uniq > femmer_matrix/$i.femmers
done
for file in *.femmers; do echo "$file"$'\n'"$(cat -- "$file")" > "$file"; done
```

## (6b) Paste FEMALES together

```bash
cd femmer_matrix
paste 15.femmers 26.femmers 34.femmers 37.femmers 32.femmers 36.femmers > Femmer_Matrix_no_MALE.txt
```


# (7) Find all femmers that are present in different numbers of females

```bash
echo "Number of Femmers in 2 or more individuals" >> Femmer_counts.txt
awk '{if ($1+$2+$3+$4+$5+$6 > 1) {print $0}}' Femmer_Matrix_no_MALE.txt | wc -l | awk '{print$1/2}' >> Femmer_counts.txt
echo "Number of Femmers in 3 or more individuals" >> Femmer_counts.txt
awk '{if ($1+$2+$3+$4+$5+$6 > 2) {print $0}}' Femmer_Matrix_no_MALE.txt | wc -l | awk '{print$1/2}' >> Femmer_counts.txt
echo "Number of Femmers in 4 or more individuals" >> Femmer_counts.txt
awk '{if ($1+$2+$3+$4+$5+$6 > 3) {print $0}}' Femmer_Matrix_no_MALE.txt | wc -l | awk '{print$1/2}' >> Femmer_counts.txt
echo "Number of Femmers in 5 or more individuals" >> Femmer_counts.txt
awk '{if ($1+$2+$3+$4+$5+$6 > 4) {print $0}}' Femmer_Matrix_no_MALE.txt | wc -l | awk '{print$1/2}' >> Femmer_counts.txt
echo "Number of Femmers in all individuals" >> Femmer_counts.txt
awk '{if ($1+$2+$3+$4+$5+$6 == 6) {print $0}}' Femmer_Matrix_no_MALE.txt | wc -l | awk '{print$1/2}' >> Femmer_counts.txt
```
