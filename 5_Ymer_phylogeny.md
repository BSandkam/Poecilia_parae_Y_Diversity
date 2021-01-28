# Building a phylogeny based on Y chromosome specific k-mers across _Poecilia parae_ males

**Goal:**
  To determine how similar Ymers are across morphs by building a phylogeny on the presence/absence of Ymers using P. picta Ymer presence/absence as a root.

General Approach:
  1) Generate file of all unique ymers (and their reverse complements) present in all parae males and in all picta
  2) Build matrix whereby we ask if each ymer is present or absent in each individual
  3) Build nexus file to run MrBayes
  4) Build phylogeny based on shared ymers
  5) Determine number of unique ymers for each clade

# 1) Generate file of all ymers from all parae individuals and picta
```bash
for Ymer_file in [files_of_kmers_not_present_in_females]
do
  cat $Ymer_file >> all_male_only_kmers
done
```


# 2a) Generate file of presence-absence of each ymer for each parae male
```bash
for Ymer_file in [files_of_kmers_not_present_in_females]
do
awk 'NR==FNR{a[$1];a[$2];next} {
  if (($1 in a) || ($2 in a))
  print "1";
else
  print "0";
}' $Ymer_file all_uniq_ymers > ymer_matrix/$Ymer_file\_ymers
done
```


## 2b) Generate file of presence-absence of each ymer for picta
```bash
awk 'NR==FNR{a[$1];a[$2];next} {
  if (($1 in a) || ($2 in a))
  print "1";
else
  print "0";
}' picta_ymers all_uniq_ymers > ymer_matrix/picta_ymers
```

## 2c) Add the name of the sample to the top of each ymer file with for ONLY ymers
```bash
cd ymer_matrix
for Ymer_presence_file in *_ymers; do echo "$Ymer_presence_file"$'\n'"$(cat -- "$Ymer_presence_file")" > "$Ymer_presence_file"; done
```

## 2d) Paste file of presence-absence of each ymer in each male together (and picta)
```bash
cd ymer_matrix
paste *_ymers > Ymer_Matrix.txt
```

## 2e) Find number of unique Parae Ymers present in at least 2 parae males
*Note: column 24 is picta and therefore not used here*
```bash
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 > 1) {print $0}}' Ymer_Matrix.txt | wc -l | awk '{print$1}' >> Number_parae_ymers.txt
```

## 2f) Generate file of presence-absence of only Ymers present in at least 2 parae males
```bash
awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 > 1) {print $0}}' Ymer_Matrix.txt > Ymer2plus_Matrix_noNames.txt
```


## 2g) Rotate to make matrix for MrBayes
*This step uses the script called matrix_rotater.awk*
```bash
awk -f matrix_rotater.awk Ymer2plus_Matrix_noNames.txt > Ymer2plus_Matrix_noNames.matrix
```

## 2h) Take taxa names from previous file and add to matrix with picta
```bash
head -1 Ymer_Matrix.txt > Ymer_Matrix_columns
awk -f matrix_rotater.awk Ymer_Matrix_columns > Ymer2plus_Matrix_Names
paste Ymer2plus_Matrix_Names Ymer2plus_Matrix_noNames.txt > Ymer2plus_Matrix_Named.matrix
```

# 3a) Make a nexus header
```bash
nano ParaeMale_2plus_wPicta.nex
#NEXUS
[ Data: Ymer2plus_Matrix_noNames.matrix ]

begin data;
  dimensions ntax=24 nchar=55901892;
  format datatype=restriction;
  matrix
```

## 3b) Add the matrix to the NEXUS header
```bash
cat Ymer2plus_Matrix_Named.matrix >> ParaeMale_2plus_wPicta.nex
```

## 3c) Make a NEXUS footer
```bash
nano end
;
end;
```

## 3d) Add the NEXUS footer to the end of the NEXUS
```bash
cat end >> ParaeMale_2plus_wPicta.nex
```


## 4) Run Mrbayes
*This runs interactively*
```bash
mb
execute ParaeMale_2plus_wPicta.nex
lset nst=1 rates=equal
mcmc ngen=50000 samplefreq=1000 printfreq=500 diagnfreq=1000
sump
sumt
```

A tree is produced that can then be read to find which samples are in which clades.


# 5a) Find the number of unique Ymers in each clade

Here we create a file called "number_of_UNIQ_ymers_for_tree.txt" and then determine how many unique Ymers are present in each clade by summing (presence=1, absence=0) across columns of our matrix (samples) and determining how many rows (Ymers) are only present in the members of the clade on the phylogeny.

*Key to sample by column in Ymer_Matrix_no_fem_wPicta.txt*

$1  $2  $3  $4  $5  $6  $7  $8  $9  $10 $11 $12 $13 $14 $15 $16 $17 $18 $19 $20 $21 $22 $23 $24
P01 P02 P03 P04 P05 P06 P09 P10 P11 P12 P13 P14 P17 P19 P20 P21 P22 P23 P24 P27 P29 P30 P31 picta

```bash
touch number_of_UNIQ_ymers_for_tree.txt
#Melanzona samples  01 02 11 12 13 14 17 19 20 21 22 23 24
echo "In all melanzona" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$2+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19 == 13) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 13)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P12 P01" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$10 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P12 P01 P02" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$10+$2 == 3) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 3)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P12 P01 P02 P19" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$10+$2+$14 == 4) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 4)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P11 P20" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($9+$15 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P12 P01 P02 P19 P11 P20" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$10+$2+$14+$9+$15 == 6) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 6)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P14 P23" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($12+$18 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P14 P23 P24" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($12+$18+$19 == 3) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 3)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P14 P23 P24 P13" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($12+$18+$19+$11 == 4) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 4)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P12 P01 P02 P19 P11 P20 P14 P23 P24 P13" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$10+$2+$14+$9+$15+$12+$18+$19+$11 == 10) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 10)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P17 P22" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($13+$17 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P17 P22 P21" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($13+$17+$16 == 3) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 3)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt

#Immaculata samples  05 06 09 10
echo "In all Immaculata" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($5+$6+$7+$8 == 4) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 4)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P09 P10" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($7+$8 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P09 P10 P06" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($7+$8+$6 == 3) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 3)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt

#Parae samples  03 04 27 29 30 31
echo "In all Parae" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($3+$4+$20+$21+$22+$23 == 6) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 6)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P29 P30" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($21+$22 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P03 P31" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($3+$23 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P04 P27" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($4+$20 == 2) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 2)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
echo "In P03 P31 P04 P27" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($3+$23+$4+$20 == 4) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 4)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt

#All Parae and Melanzona
echo "In all Parae and Melanzona" >> number_of_UNIQ_ymers_for_tree.txt
awk '{if ($1+$2+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$3+$4+$20+$21+$22+$23 == 19) {print $0}}' Ymer_Matrix.txt | awk '($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 == 19)' | wc -l | awk '{print$1}' >> number_of_UNIQ_ymers_for_tree.txt
```


## 5b) Find number of shared P. parae Ymers that are also Ymers in P. picta

Here we create a file called "number_of_Picta_Ymers_by_sample.txt" and then determine how many Ymers are present in picta and each male.

```bash
awk '{if ($24 == 1) {print $0}}' Ymer_Matrix.txt | awk '{if ($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24 > 2) {print $0}}' > picta_w_parae2_matrix
echo "In Picta and P01" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($1+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P02" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($2+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P03" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($3+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P04" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($4+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P05" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($5+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P06" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($6+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P09" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($7+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P10" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($8+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P11" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($9+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P12" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($10+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P13" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($11+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P14" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($12+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P17" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($13+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P19" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($14+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P20" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($15+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P21" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($16+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P22" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($17+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P23" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($18+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P24" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($19+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P27" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($20+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P29" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($21+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P30" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($22+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
echo "In Picta and P31" >> number_of_Picta_Ymers_by_sample.txt
awk '{if ($23+$24 == 2) {print $0}}' picta_w_parae2_matrix | wc -l | awk '{print$1}' >> number_of_Picta_Ymers_by_sample.txt
```
