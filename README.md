# Primer Design tool

Purpose of this tool is to design (nested) primers for multiplex PCRs and RTs. 

Input can be a list of gene names or a list of genomic coordinates to target.

Primers are checked against each other to select sets of primers with minimal 3' complementary so as to avoid formation of dimers. Primers are further blasted against the genome or transcriptome (depending on use case) to make sure no unintended products form. 

Several use cases are described in detail below.

## Basic setup

Make sure that the following R libraries are avilable:

```
library(getopt)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(magrittr)
library(stringr)
```
and, depending on the genome version you wish to use:

```
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
```
Furthermore, I currently have a number of resources hardcoded into the Script, these should be bundled upon release.

Genome and transcriptome FASTA files (for construction of blast DB)

```
/g/steinmetz/genome/Homo_sapiens/37.68/fasta/DNA/Homo_sapiens.GRCh37.68.dna.chromosomes.withERCC.fa
/g/steinmetz/genome/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.ERCC.fa
/g/steinmetz/genome/Homo_sapiens/38.79_tx/fasta/Human_ERCC_combined.convenient.fa

/g/steinmetz/genome/Mus_musculus/38.73/fasta/Mus_musculus.GRCm38.73.dna.chromosome.all.spikes.fa
/g/steinmetz/genome/Mus_musculus/38.91/Mus_musculus.GRCm38.cdna.convenient.all.fa
```

GFFs of protein coding exons, for construction of a TxDb object

```
/g/steinmetz/project/singcellTxn/CRISPRdrop/hg38_annotation/Homo_sapiens.GRCh38.89.chr.protein-coding.exon.gtf.gz
/g/steinmetz/velten/Scripts/Primer/annotation/Homo_sapiens.GRCh37.87.chr.gtf.gz
/g/steinmetz/velten/Scripts/Primer/annotation/Mus_musculus.GRCm38.92.chr.gtf.gz
```

Executables

```
/g/steinmetz/gschwind/software/blast+/ncbi-blast-2.6.0+/bin
/g/software/linux/pack/primer3-2.3.4/bin/primer3_core
/g/steinmetz/velten/Scripts/Primer/parsePrimers.pl
```
## Use cases

### Use case 1a: Target a list of genes during reverse transcription

Example: You have a list of marker genes that you need to confidently identify in each cell in a smart-seq2 experiment. The input gene list should be provided as a txt file with one gene symbol per line. 
Mouse gene symbols are in the format `Actb`
Human gene symbols are in the format `ACTB`

Primers will always be place on exons present in all refseq coding transcripts, if the RT product is then still longer than `--lengthRTProduct`. 
Often, no exon is present in all refseq transcripts. The program will then add a warning to the final CSV output and use the exons it considers best.

Usage example (see below for all options)

```
Rscript Masterscript.R --geneSymbols Markergenes.txt --cDNA --genome mm10 \
--nested RTonly --lengthRTProduct 300 --csv targeted_rt_primers.csv --bed targeted_rt_primers.bed --blast --optim --prefixOuter AAGCAGTGGTATCAACGCAGAGT
#RT primers are considered outer primers in the context of all arguments - e.g. --prefixOuter preifxes all RT primers with a common sequence (here, the ISPCR sequence)
```

### Use case 2: Target genomic coordinates during reverse transcription

Example: You are looking at some cancer mutations and want to increase their coverage in a smart-seq2 experiment. The input should be a list of genomic coordinate provided as a txt file, in format `chr2:25470536`
Primers will always be place on the same exon carrying the mutation or on another exon present in all transcripts containing the exon of interest (let's double check this behavior!)

```
Rscript Masterscript.R --coordinates coordinates.txt --cDNA --genome hg19 \
--nested RTonly --lengthRTProduct 600 --csv targeted_rt_primers.csv --bed targeted_rt_primers.bed --blast --optim --prefixOuter AAGCAGTGGTATCAACGCAGAGT
```

### Use case 3: Target genomic coordinates during reverse transcription + smartSeq PCR

As in 10.1038/nm.4336

```
Rscript Masterscript.R --coordinates coordinates.txt --cDNA --genome hg19 \
--nested RT --lengthRTProduct 1000 --lengthInnerProduct 500,1000 --csv targeted_rt_primers.csv --bed targeted_rt_primers.bed --blast --optim --prefixInner AAGCAGTGGTATCAACGCAGAGT
#RT primers are considered outer primers and PCR primers are considered inner primers in that setting

```


### Use case 4: Target genomic coordinates for nested PCRs

Genomic coordinates can be on gDNA (e.g. when amplifying material from single cells or single-cell derived colonies) or cDNA (e.g. when amplifying from a smartseq2 library).

```
Rscript Masterscript.R --coordinates coordinates.txt --gDNA --genome hg19 \
--nested full --csv targeted_nested_primers.csv --bed targeted_nested_primers.bed --blast --optim --prefixInner XXXXX
#RT primers are considered outer primers and PCR primers are considered inner primers in that setting

```

## Overview of all options

| Flag        | Argument       | Description  |
| ------------- |:-------------:| -----:|
|  --geneSymbols  | file | file containing Gene Symbols |
|  --coordinates    | file      |   file containing genomic coordinates, mutually exclusive with geneSymbols |
| --nested | RT, RTonly, none or full | Whether to select PCR primers (none), nested PCR primers (full), an outer RT primer + inner PCR primers (RT) or an RT primer only. |
| --csv | file | csv file for output |
| --bed | file | bed file for output (optional). Creates a bed file highlighting the position of the primers finally selected as well as all off-target sites they bind to. |
| --gDNA | none  |  Input is gDNA (not compatible with gene symbols) |
| --cDNA | none  |  Input is cDNA |
| --genome | hg19, hg38 or mm10 | Genome version to use |
| --optim | none | Select set of primer pairs with minimal 3' overlap to reduce byproduct formation. This slows down computations a lot if --numprimers or the number of targets is too high |
| --blast | none | Blast all primer candidates to select candidates that do not form byproduct from. Note that final primers will always be blasted aganst the genome if bed output is desired |
| --genomedb | none | Path to Blast database for genome. If left out, will be created if needed|
| --txdb | none | Path to Blast database for transcriptome. If left out, will be created if needed|
| --lengthInnerProduct | int,int | Length of inner PCR product (min, max) |
| --lengthOuterProduct | int,int | Length of outer PCR product (min, max) |
| --lengthRTProduct | int | Minimum length of RT product (if possible) |
| --TmOuter | num,num,num | Melting temperature of outer or RT primers (optimal,min,max) |
| --TmInner | num,num,num | Melting temperature of inner primers (optimal,min,max) |
| --nCharOuter | int,int,int |Length of outer or RT primers (optimal,min,max) |
| --nCharInner | int,int,int | Length of inner primers (optimal,min,max) |
| --prefixOuter | DNA sequence | Prepend all outer or RT primers with a DNA sequence |
| --prefixInner | DNA sequence | Prepend all inner primers with a DNA sequence |
| --numprimers | int | Number of primers to spit out during initial primere3 calls - higher number make blasting and optimisation slow but increase the chance to find good primers |
| --help | none | Prints a list of arguments |
| --verbose | none | Create a lot of output |
