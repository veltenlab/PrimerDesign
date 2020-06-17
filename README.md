# Primer Design tool

Purpose of this tool is to design (nested) primers for multiplex PCRs and RTs. 

Input can be a list of gene names or a list of genomic coordinates to target.

Primers are checked against each other to select sets of primers with minimal 3' complementary so as to avoid formation of dimers. Primers are further blasted against the genome or transcriptome (depending on use case) to make sure no unintended products form. 

Several use cases are described in detail below.

Run time depends on the number of primers tested but typically stays below 30 minutes on a default laptop computer.

## Requirements and first run

PrimerDesign tool runs under Linux and requires [primer3](https://github.com/primer3-org/primer3), [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), and [cliquer](https://users.aalto.fi/~pat/cliquer.html) to be available on your system; the paths need to be specified by the user in the first run of the program. All other requirements are installed automatically when the package if first run; in case you do not have write permissions to your `.libPaths()`, uncomment the first line of the scripot and specifiy a writeable folder that R libraries will be installed to.

Genome data will be downloaded and prepared in the first run for a given genom (mm10, hg19 or hg38), so please be patient when you first run this software. In total, installation may take between 5 minutes and an hour depending on the prior availability of packages and internet connection.

## Use cases


### Use case 1: Target genomic coordinates for MutaSeq

Example: You are looking at some cancer mutations and want to increase coverage. The input should be a list of genomic coordinate provided as a txt file, in format chr2:25470536.

```
Rscript Masterscript.R --coordinates coordinates.txt --cDNA --genome hg19 \
--readlength 75  --nested none --csv targeted_MutaSeq_primers.csv --bed targeted_MutaSeq_primers.bed --blast --optim

```

### Use case 2: Target genomic coordinates for nested PCR amplification

Genomic coordinates can be on gDNA (e.g. when amplifying material from genomic DNA) or cDNA (e.g. when amplifying from a cDNA library).

```
Rscript Masterscript.R --coordinates coordinates.txt --gDNA --genome hg19 \
--readlength 75 --nested full --csv targeted_nested_primers.csv --bed targeted_nested_primers.bed --blast --optim

```

### Use case 3: Target genes during reverse transcription

Example: You have a list of marker genes that you need to confidently identify in each cell in a smart-seq2 experiment. The input gene list should be provided as a txt file with one gene symbol per line. 
Mouse gene symbols are in the format `Actb`
Human gene symbols are in the format `ACTB`

Primers will always be place on exons present in all refseq coding transcripts, if the RT product is then still longer than `--lengthRTProduct`. 
Often, no exon is present in all refseq transcripts. The program will then add a warning to the final CSV output and use the exons it considers best.

Usage example (see below for all options)

```
Rscript Masterscript.R --geneSymbols Markergenes.txt --cDNA --genome mm10 \
--nested RTonly --lengthRTProduct 300 --csv targeted_rt_primers.csv --bed targeted_rt_primers.bed --blast --optim 
```

Alternatively, --nested RT will additionally design PCR primers.


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
| --optim | none | Select set of primer pairs with minimal 3' overlap to reduce byproduct formation.|
| --optimThr | int | Threshold of primer complimentarity to allow. Defaults to 15. |
| --optimThrSame | int | Threshold of self-complimentarity to allow. Defaults to 30. |
| --blast | none | Blast all primer candidates to select candidates that do not form byproduct from. Note that final primers will always be blasted aganst the genome if bed output is desired |
| --genomedb | none | Path to Blast database for genome. If left out, will be created if needed|
| --txdb | none | Path to Blast database for transcriptome. If left out, will be created if needed|
| --readlength | int | Length of the first sequencing read to ensure mutation is covered |
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
