**Goals**
  To determine relative coverage of males and females when aligned to the best female de novo assembly (P32) which was used in the RACA pipeline

Approach:
  1) Align all trimmed Illumina reads (both 10X pipeline and straight Illumina) from all samples to P32 female
  2) Determine scaffold coverage differences between each male morph and females
  3) Determine position and coverage differences of the scaffolds that RACA assembled into chromosomes so they can be graphed


# Step 1) Align all raw reads (both 10X pipeline and straight Illumina) from all samples to the female P32 de novo scaffolds greater than 10kb (used for the RACA pipeline)

## 1a) Index the P32_greater_10kb file
Parameters
  bwa
    index = indexes the de novo genome for alignment

```bash
bwa index P32_greater_10kb.fa
```


## 1b) Align fastq reads from true illumina to P32_greater_10kb
Parameters
  bwa aln
    -t = 20 threads

*Run for each sample*
```bash
bwa aln -t 20 P32_greater_10kb.fa [forward_paired_reads.fastq] > [sample_for.sai]
bwa aln -t 20 P32_greater_10kb.fa [reverse_paired_reads.fastq] > [sample_rev.sai]
```


## 1c) Generate SAM files from the .sai files
Parameters
  bwa sampe

*Run for each sample*
```bash
bwa sampe P32_greater_10kb.fa [sample_for.sai] [sample_rev.sai] [forward_paired_reads.fastq] [reverse_paired_reads.fastq] > [sample.sam]
```


## 1d) Find only uniquely mapping reads for further analysis

*Run for each sample*
```bash
grep 'XT:A:U' [sample.sam] > [sample.uniq.sam]
```



# Step 2) Determine scaffold coverage differences
## 2a) Estimate coverage with [soap.coverage v2.7.7](http://soap.genomics.org.cn)
Parameters
  soap.coverage
    -cvg = sequencing coverage mode
    -refsingle = input reference fasta file
    -i = input filename
    -sam = input is in sam format
    -o = output filename
    -plot [filename] 0 1000 = plots coverage from 0 to 1000
    -p = number of threads

*Run for each sample*
```bash
soap.coverage -cvg -refsingle P32_greater_10kb.fa -i [sample.uniq.sam] -sam -o [sample_soapCov.txt] -plot [sample_distribution.txt] 0 1000 -p 24
```


## 2b) Prepare soap.coverage output files for downstream analysis
- Remove the "Percentage:" and the "Depth:" from each file.

```bash
sed -i 's/Depth://g' *soapCov.txt
sed -i 's/Percentage://g' *soapCov.txt
```


## 2c) Move male coverage files of each sample to respective morph and total male directories, move female coverage files to female directory
```bash
cp [Red_sample_soapCov.txt] red/
cp [Yellow_sample_soapCov.txt] yellow/
cp [Blue_sample_soapCov.txt] blue/
cp [Immaculata_sample_soapCov.txt] immac/
cp [Parae_sample_soapCov.txt] parae/
cp [Male_sample_soapCov.txt] male/
cp [Female_sample_soapCov.txt] female/
```


## 2d) Use the extract_coverage.py script to determine coverage of each morph

```bash
python extract_coverage.py red/ red_coverage.txt
python extract_coverage.py yellow/ yellow_coverage.txt
python extract_coverage.py blue/ blue_coverage.txt
python extract_coverage.py immac/ immaculata_coverage.txt
python extract_coverage.py parae/ parea_coverage.txt
python extract_coverage.py male/ males_coverage.txt
python extract_coverage.py female/ females_coverage.txt
```


## 2e) Use the foldchange_coverage.py script to calculate M:F log average coverage for each morph and for all males

```bash
python foldchange_coverage.py females_coverage.txt red_coverage.txt redF_coverage_fold_change.txt
python foldchange_coverage.py females_coverage.txt yellow_coverage.txt yellowF_coverage_fold_change.txt
python foldchange_coverage.py females_coverage.txt blue_coverage.txt blueF_coverage_fold_change.txt
python foldchange_coverage.py females_coverage.txt parae_coverage.txt paraeF_coverage_fold_change.txt
python foldchange_coverage.py females_coverage.txt immaculata_coverage.txt immaculataF_coverage_fold_change.txt
python foldchange_coverage.py females_coverage.txt males_coverage.txt maleF_coverage_fold_change.txt
```



# 3) Determine position and coverage of scaffolds that RACA assembled into chromosomes
## 3a) Find RACA placed scaffolds with counter_RACAscaff.py script
*Note: need the rec_chrs.ppar.segments.refined.txt file of parae segments from the RACA output in the Out_RACA directory*

```bash
python counter_RACAscaff.py [rec_chrs.ppar.segments.refined.txt] ppar RACA_segments_used.txt
```


## 3b) Find where the RACA placed scaffolds are located relative to the Xiph genome using scaff_posinfo_RACA.py script
*Note: need the Conserved.Segments file from the RACA output in the Out_RACA/SFs directory*
```bash
python scaff_posinfo_RACA.py ppar Conserved.Segments RACA_segments_used.txt [rec_chrs.ppar.segments.refined.txt] RACA_scaff_posinfo.txt
```


## 3c) Make file of uniquely mapping scaffolds and their coverage differences between each morph and females with extract_cov_raca_scaffs.py script

```bash
python extract_cov_raca_scaffs.py RACA_scaff_posinfo.txt redF_coverage_fold_change.txt red_scaf_cov.txt
python extract_cov_raca_scaffs.py RACA_scaff_posinfo.txt yellowF_coverage_fold_change.txt yellow_scaf_cov.txt
python extract_cov_raca_scaffs.py RACA_scaff_posinfo.txt blueF_coverage_fold_change.txt blue_scaf_cov.txt
python extract_cov_raca_scaffs.py RACA_scaff_posinfo.txt paraeF_coverage_fold_change.txt parae_scaf_cov.txt
python extract_cov_raca_scaffs.py RACA_scaff_posinfo.txt immaculataF_coverage_fold_change.txt immaculata_scaf_cov.txt
python extract_cov_raca_scaffs.py RACA_scaff_posinfo.txt maleF_coverage_fold_change.txt male_scaf_cov.txt
```
