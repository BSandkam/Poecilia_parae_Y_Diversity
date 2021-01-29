**Goal**
  To compare percent of sequence comprised of interspersed repeats between Y linked scaffolds of each morph and the rest of the genome (indicates transposible elements (TE) presence and activity)

Approach
  1) Generate database of repeat sequence from P. parae sequence
    a) Run RepeatModeler on female de novo genome to generate models for X linked and autosomal sequence repeats
    b) Combine the female Repeat database with the repeat databases built from the morph specific Y-linked scaffolds from males
  2) Run RepeatMasker using the combined repeat database on females and the morph specific Y-linked scaffolds


# 1) Generate database of repeat sequence from P. parae sequence

## 1a) Run RepeatModeler on female de novo genome to generate models for X linked and autosomal sequence repeats

*The next step uses the python script fastShortHeads.py to shorten the fasta headers/id to work with RepeatModeler*

Parameters
  seqtk seq
    -l 60 = writes 60 bases per line
    -L 5000 = drops sequences less than 5000 bp (causes errors for RepeatModeler - full sequences are used for RepeatMasking later in the pipeline)
  fasShortHead.py
    -l 20 = number of characters to shorten fasta header names to
  BuildDatabase (part of the RepeatModeler package)
    -engine ncbi = search engine to use
    -name makeranno = name of database to create
  RepeatModeler
    -rand 12345 = setting the random number generator so runs are repeatable
    -engine ncbi = search engine to use
    -pa 34 = 34 threads
    -database female = name of database to use

```bash
seqtk seq -l 60 -L 5000 P32_greater_10kb.fa | python fasShortHead.py -l 20 - > genome2modeler.fas

BuildDatabase -engine ncbi -name female P32_greater_10kb.fa

RepeatModeler -srand 12345 -engine ncbi -pa 34 -database female
```
*Note: this creates as -families.fa file to be used in next step*


## 1b) Run RepeatModeler on morph-linked scaffolds from each male

Parameters
  seqtk seq
    -l 60 = writes 60 bases per line
    -L 5000 = drops sequences less than 5000 bp (causes errors for RepeatModeler - full sequences are used for RepeatMasking later in the pipeline)
  fasShortHead.py
    -l 20 = number of characters to shorten fasta header names to
  BuildDatabase (part of the RepeatModeler package)
    -engine ncbi = search engine to use
    -name makeranno = name of database to create
  RepeatModeler
    -rand 12345 = setting the random number generator so runs are repeatable
    -engine ncbi = search engine to use
    -pa 34 = 34 threads
    -database female = name of database to use

```bash
for i in `ls *.out | sed 's/_mer_scaffs_5plus_all.out//g'`
do
seqtk seq -l 60 -L 5000 $i\_mer_scaffs_5plus_all.out > $i\_for_modeling.fa
BuildDatabase -engine ncbi -name $i $i\_for_modeling.fa
RepeatModeler -srand 12345 -engine ncbi -pa 34 -database $i
done

cp *-families.fa ../repeat_families/
```


## 1c) Combine the female Repeat database with the repeat databases built from male morph specific Y-linked scaffolds and the full Actinopterygii database
*Note: Before proceeding move each of the -families.fa files created by running RepeatModeler on the male morph-linked Y scaffolds to the same directory as the female*
*The queryRepeatDatabase.pl script is part of the RepeatMasker utility package*
```bash
cat *-families.fa | awk '{print $1}' > _tmp.families.fa
queryRepeatDatabase.pl -species Actinopterygii > Actinopterygii.rmLib.fa
cat _tmp.families.fa Actinopterygii.rmLib.fa > combined.RM.lib
```


# 2) Run RepeatMasker using the combined repeat database on each of the females and the morph specific Y-linked scaffolds from each male
*Note: move the scaffold files to the same folder, and add the file of the combined repeat masker library*

In total run individually on:
  Each of the 6 Female de novo genomes
  All scaffolds containing more than 5 of the respective morph-mers (morph specific Y-linked kmers) for each of the 23 male individuals
    - 13 melanzona
    - 6 parae morph
    - 4 immaculata

*Note: comparing morph specific scaffolds to female genomes because we were conservative with morph-specific scaffold identification and therefore may not have identified all the morph unique sequence. By comparing to the female genome we are sure we are comparing morph-specific Y linked sequence to autosomal and X linked sequence.*

Parameters
  RepeatMasker
    -engine ncbi = search engine to use
    -pa 30 = 30 threads
    -noisy = print status to stdout while running
    -xsmall = return repetitive regions as small characters rather than masked
    -alignments = write alignments in .align output files
    -dir = writes output to this directory
    -lib combined.RM.lib = use the custom library created from P. parae female genome, morph-linked Y scaffolds, and the full Actinopterygii database

```bash
for i in *.fas
do
  RepeatMasker -engine ncbi -pa 30 -noisy -xsmall -alignments -dir . -lib combined.RM.lib $i
done
```
