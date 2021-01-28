# Pipeline to Annotate Y-linked contigs

## Pre-requisites

- [NextFlow](https://www.nextflow.io) (tested with version 19.04.1)
- [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) (`RepeatModeler`, `BuildDatabase`) (tested with version open-1.0.10)
- [RepeatMasker](http://repeatmasker.org) (`RepeatMasker`, `queryRepeatDatabase.pl`, `rmblast`) (tested with version open-4.0.7)
- [BUSCO](https://busco.ezlab.org) (the executable should be called `busco`) (tested with version 3.0.2)
    + with the necessary OrthoDB database (for example `Actinopterygii`)
- [SNAP](https://github.com/KorfLab/SNAP) (`snap` and all associated scripts) (tested with version 2006-07-28)
- [Augustus](http://bioinf.uni-greifswald.de/augustus/) (`augustus`) (tested with version 3.2.3)
    + the `$AUGUSTUS_CONFIG_PATH` variable must be set and the path writeable without root permissions
- [MAKER](https://www.yandell-lab.org/software/maker.html) (`maker` and all associated scripts) (tested with version 2.31.10)
- [BEDtools](https://bedtools.readthedocs.io/en/latest/) (`bedtools`)
- [seqtk](https://github.com/lh3/seqtk)
- Python3 (tested with Python 3.6, but should work with other 3.X versions) with the following non-standard libraries


## Datasets

1. fastA file with contigs to annotate
    - fastA headers should not contain weird characters like `|`. For example, if you're running this from the representative sequences identified with CD-HIT and the fasta headers have `|` in their string, replace it with other character such as `_`.
2. proteins
    - full protein amino acid sequences, in fastA format, for 8 other fish species downloaded from Ensembl (see table below).
3. EST transcripts in fastA format
    - EST data was obtained from Dreyer et al. BMC Genomics 2007, and used 3421 ESTs matching the word "testis". Fasta files for all EST were concatenated into a single file `Dreyer_etal_2007.seqs.testis.fas`
4. reference-guided transcriptome assembly created with HISAT and StringTie (but other tools should also work)
    - two RNA-seq datasets from Sharma et al. 2014, both from male samples, one from testis and another from embryo. However, this was done using the guppy reference genome, so might not be appropriate for other species.

Protein datasets from 8 fish species:

| Species                | Common Name              |
|:-----------------------|:-------------------------|
| Danio rerio            | Zebrafish                |
| Gasterosteus aculeatus | three-spined stickleback |
| Oryzias latipes Hd-rR  | Medaka Hd-rR             |
| Poecilia latipinna     | Sailfin molly            |
| Poecilia mexicana      | Shortfin molly           |
| Poecilia reticulata    | Guppy                    |
| Takifugu rubripes      | Fugu                     |
| Xiphophorus maculatus  | Southern platyfish       |


## Preparation

Assuming all **protein** fastA files reside in a folder named `Proteins`, create a file of filenames `proteins.fof` with the relative path to the fasta files, one per line. For example:

```
../Proteins/Danio_rerio.GRCz11.pep.all.fa
../Proteins/Gasterosteus_aculeatus.BROADS1.pep.all.fa
../Proteins/Oryzias_latipes.ASM223467v1.pep.all.fa
../Proteins/Poecilia_latipinna.P_latipinna-1.0.pep.all.fa
../Proteins/Poecilia_mexicana.P_mexicana-1.0.pep.all.fa
../Proteins/Poecilia_reticulata.Guppy_female_1.0_MT.pep.all.fa
../Proteins/Takifugu_rubripes.FUGU5.pep.all.fa
../Proteins/Xiphophorus_maculatus.X_maculatus-5.0-male.pep.all.fa
```
*Note: Takifugu_rubripes.FUGU5.pep.all.fa was taken from http://www.magic.re.kr/publicdb/takifugu_rubripes/pep/ instead*

Similarly, create second file of filenames `est.fof` with the relative path to the fasta files with EST and transcripts, for example:

```
../EST/Dreyer_etal_2007.seqs.testis.fas
../EST/male_embryo.libmerge.transcripts.fpkm1.ncrna.fna
../EST/male_testis.libmerge.transcripts.fpkm1.ncrna.fna
```


## Make a copy of the fasta file of morph specific scaffolds that has header with no pipes
*Note: Run for each morph*
```bash
sed 's/\|/_/g' [morph_scaffolds.fas] > [morph_scaffolds_noPipe.fas]
```

## Run the pipeline
```bash
nextflow run makerAnno.nf \
    --fas [morph_scaffolds_noPipe.fas] \
    --est est.fof --prot proteins.fof \
    --rm_mod --seqL 1000 \
    --rm_taxa Actinopterygii \
    --line ../actinopterygii_odb9\
    --spc zebrafish \
    --cpu 30 --prefix [morph_mers]
```

- `--fas` takes the fastA file with the contigs to annotate
- `--est` and `--prot` tale the file of filenames created above for EST/transcripts and proteins
- `--rm_mod` tells the pipeline to run RepeatModeler using only contigs at least 1kb long (`--seqL 1000`). This filter only applies to RepeatModeler and not to RepeatMasker or the other steps
- `--rm_taxa Actinopterygii` tells RepeatMasker to use the actinopterygii repeat database from RepBase (note that RepBase should have been installed with RepeatMasker)
- `--line ~/DB/BUSCO/actinopterygii_odb9` tells BUSCO to use the actinopterygii OrthoDB database
- `--spc zebrafish` tells Augustus to use the HMM training sets of zebrafish

Everything else in handled automatically by the pipeline. `makerAnno.nf` is simply a text file so it's possible to open it a text editor and see all the steps.


## Post-processing

When the pipeline finishes, there are several post-processing steps, the first of which is to organize all output files from MAKER. For this install the tools from `https://github.com/NBISweden/GAAS` following the instructions to add them to both `PERL5LIB` and `PATH` environmental variables.

Make the following directories
```bash
cd ../Annotation/running
mkdir 05.MAKER_OUTPUTS
mkdir 05.MAKER_OUTPUTS/01.BLAST_Uniprot
mkdir 05.MAKER_OUTPUTS/02.InterProScan
```

# Map IDs from the maker output
*This will let us rename the maker generated sequences in the fasta and gff files in the next step*

```bash
maker_map_ids --prefix [morph_mers_] --justify 8 [morph_mers.round2.all.maker.noseq.gff] > [morph_mers.round2.maker.rename.map]
```

# Rename sequences to be shorter than maker output names

```bash
#Rename sequences in fasta
for f in *.fasta; do
	cp $f ${f/.fasta/.rename.fasta}
	map_fasta_ids [morph_mers.round2.maker.rename.map] ${f/.fasta/.rename.fasta}
done

#Rename sequences in gff
for g in *.gff; do
	cp $g ${g/.gff/.rename.gff}
	map_gff_ids [morph_mers.round2.maker.rename.map] ${g/.gff/.rename.gff}
done

#Rename the gff file according to new names
map_gff_ids [morph_mers.round2.maker.rename.map] [morph_mers.round2.all.maker.rename.gff]
```

# Prepare the uniprot database
First Download the uniprot database ( https://www.uniprot.org/downloads - downloaded the "Reviewed (Swiss-Prot)" in fasta format on 2020-04-01)

Index the uniprot database
```bash
index uniprot db
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```

# BLAST to uniprot database

```bash
blastp -task blastp -query [morph_mers.round2.all.maker.proteins.rename.fasta] -db uniprot_sprot.fasta -out maker_annotation.proteins.rename_uprot.blast.out \
	-evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -num_threads 30
```

# Install iprscan5
- Downloaded from https://github.com/ebi-wp/webservice-clients and installed on 2020-04-01

# Run iprscan5
```bash
#Start virtual environment
virtualenv -p `which python` env
source ../webservice-clients-master/env/bin/activate

#Run iprscan5 - change email
iprscan5.py --sequence [morph_mers.round2.all.maker.proteins.rename.fasta] --multifasta --useSeqId \
	--email '_@_._' --maxJobs 10 \
	--title Parae-ymers --outformat tsv --appl PfamA

echo -n > [morph_mers.round2.all.maker.proteins.rename_iproscan.out]

#Deactivate virtual environment
deactivate
```

# Clean up files
Combine iprscan5 output files
```bash
for f in *.tsv.txt; do cat $f >> [morph_mers.round2.all.maker.proteins.rename_iproscan.out] ; done
```

Run maker_functional_fasta
```bash
maker_functional_fasta uniprot_sprot.fasta maker_annotation.proteins.rename_uprot.blast.out [morph_mers.round2.all.maker.proteins.rename.fasta] > [morph_mers_maker_annotation.proteins.rename.funcAnno.faa]
```

Run maker_functional_gff
```bash
maker_functional_gff uniprot_sprot.fasta maker_annotation.proteins.rename_uprot.blast.out [morph_mers.round2.all.maker.rename.gff] > [morph_mers_maker_annotation.rename.funcAnno.gff]
```

Run ipr_update_gff
```bash
ipr_update_gff [morph_mers_maker_annotation.rename.funcAnno.gff] [morph_mers.round2.all.maker.proteins.rename_iproscan.out] > [morph_mers.maker_annotation.rename.funcAnno.iproscan.gff]
```


Grab gene names from gff file
```bash
sed 's/\%/_/g' [morph_mers.maker_annotation.rename.funcAnno.iproscan.gff] > [morph_gene_names.gff]
awk -F"\t" '$3=="gene" {for (x=nr; x<NF; x++) {printf $x " "} print $NF }' nr=11 [morph_gene_names.gff} | awk -F" " '{for (x=nr; x<NF; x++) {printf $x " "} print $NF }' nr=3
```
