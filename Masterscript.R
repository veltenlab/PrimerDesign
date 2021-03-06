## set your R library path source here - defaults to .libPaths()
###.libPaths("/g/steinmetz/velten/Software/RPacks3.4.0/")

options(warn=-1)
if (!interactive()){
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename,"check.R"))

## create source directory for genome files
if(!dir.exists(file.path(script.basename,"genomeFiles"))){
    cat("Creating directory for genomes 'genomesFiles'\n")
    dir.create(file.path(script.basename,"genomeFiles"))
}

###### 1. User input ######

spec = matrix(c(
  'geneSymbols', 's', 1, "character",
  'coordinates' , 'x', 1, "character",
  'gDNA' , 'g', 0, "logical",
  'cDNA', 'c', 0, "logical",
  'optim', 'y',0,"logical",
  'optimThr','q',1,"integer",
  'optimThrSame','Q',1,"integer",
  'optimOuter','Y',0,"logical",
  'blast', 'z', 0, "logical",
  'genomedb', 'E',1,'character',
  'txdb', 'F',1,'character',
  'genome' ,'G', 1, "character",
  'nested' , 'n', 1, "character",
  'lengthInnerProduct', 'I',1,"character",
  'lengthOuterProduct', 'O',1,"character",
  'lengthRTProduct', 'R',1,"integer",
  'TmOuter',"T",1,"character",
  'TmInner',"S",1,"character",
  'nCharOuter',"B",1,"character",
  'nCharInner',"D",1,"character",#  'prefixInner',"p",1,"character",
  'prefixInnerLeft',"l",1,"character",
  'prefixInnerRight',"r",1,"character",
  "prefixOuter","P",1,"character",
  'numprimers',"N",1,"integer",
  'readlength','L',1,"integer",
  'csv','o',1,"character",
  'bed','b',1,"character",
  'help','h',0,"logical",
  'verbose','v',0,'logical'
), byrow=TRUE, ncol=4)
 opt = getopt(spec)
 save(opt , script.name, script.basename, file = "opt.rda")
 
if ( !is.null(opt$help) ) {
   cat(getopt(spec, usage=TRUE))
  q(status=1)
}
 } else {
   load("opt.rda")
   cat("Restoring previous parameters\n")
   source(file.path(script.basename,"check.R"))
 }
 
 


 
#set some reasonable defaults for the options that are needed,
#but were not specified.

if (is.null(opt$csv)) stop("Please specify an output csv file through --csv")

if (is.null(opt$verbose)) verbose <- F else verbose <- T

 #We need the option to do nested designs, either both-sided or just an RT primer in the outer set
 if (!is.null(opt$nested)) {
   nested <- opt$nested #other values: none, RT, RTonly
   if (!nested %in% c("RT","RTonly","none","full")) stop("Invalid value for --nested. Permitted values are RT, RTonly, none, full.\n")
 } else {
   nested <- "none"
   cat("Assuming that no nested design is needed...\n")
 }
 
if ( !is.null(opt$geneSymbols )  ) { 
  
  targets <- readLines( opt$geneSymbols )
  targets <- targets[!grepl("^#", targets)]
  targets <- unique(targets)
  
  input <- "cDNA"
  targetType <- "geneName"
  
  if (!is.null(opt$gDNA ) ) {
    cat("Using gene symbols, ignoring gDNA flag - using cDNA instead\n")
  }
  if (!is.null(opt$coordinates ) ) {
    cat("Gene symbols and coordinates cannot both be set - ignoring coordinates\n")
  }
} else if (!is.null(opt$coordinates )) {
  targets <- readLines( opt$coordinates )
  targets <- targets[!grepl("^#", targets)]
  targets <- unique(targets)
  targetType <- "coordinate"
  if (!is.null(opt$cDNA )) {
    input <- "cDNA"
  } else if (!is.null(opt$gDNA)) {
    input <- "gDNA"
  } else {
    input <- "gDNA"
    cat("No input type specified - using genomic DNA\n")
  }
  
} else {
  stop("You need to either input a list of gene symbols through --geneSymbols or a list of genomic coordinates through --coordinates\n")
}


if (targetType == "coordinate") {
  if (!all(grepl("^[^:]+:\\d+$",targets))) stop("Wrong input format for target type coordinate")
}

if (!is.null(opt$genome)) {
  genome <- opt$genome
} else {
  stop("You need to specify if your input is from mm10, hg19 or hg38 via the --genome flag\n")
}

margin <- 2 #when pulling out sequences as template for primer3 - by what factor should sequence be larger than target?

if (is.null(opt$lengthInnerProduct)) {
  length_inner <- c(90,145)
  cat("Inner product length not specified, setting to 90-145\n")
} else {
  length_inner <- as.integer(strsplit(opt$lengthInnerProduct,",")[[1]])
}

if (is.null(opt$lengthOuterProduct) & !nested %in% c("RT","RTonly")) {
  length_outer <- c(200,350)
  cat("Outer product length not specified, setting to 200-350\n")
} else if (!nested %in% c("RT","RTonly")) {
  length_outer <- as.integer(strsplit(opt$lengthOuterProduct,",")[[1]])
}

if (!is.null(opt$lengthRTProduct)) {
  if (nested %in% c("RT","RTonly")) {
    length_rt <- opt$lengthRTProduct
    length_outer <- rep(length_rt/margin,2)
  } else {
    cat("Nested is not RT or RTonly - ignoring --lengthRTProduct\n")
  }
} else if (nested %in% c("RT","RTonly")) {
  length_rt <- margin * max(length_inner)
  length_outer <- rep(length_rt/margin,2)
}

if(is.null(opt$prefixInnerLeft)) prefixInnerLeft <- "" else prefixInnerLeft <- opt$prefixInnerLeft
if(is.null(opt$prefixInnerRight)) prefixInnerRight <- "" else prefixInnerRight <- opt$prefixInnerRight
if (prefixInnerLeft != "" | prefixInnerRight != "") cat("Prepending", prefixInnerLeft, "to the primer that is CLOSER to the mutation\n")



source(file.path(script.basename,"accessory_functions_TxDb.R"))
source(file.path(script.basename,"blast_primers.R"))

suppressMessages({
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(magrittr)
library(stringr)
})

cat("Preparing Genomes...\n")

if (genome == "mm10") {
  suppressMessages({
    pkgTest("BSgenome.Mmusculus.UCSC.mm10", bc=T)
    pkgTest("org.Mm.eg.db", bc=T)
  })
  Bsgenome <- BSgenome.Mmusculus.UCSC.mm10
  
  if(!file.exists(file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.92.chr.gtf.gz"))) {
    cat("This appears to be the first time that you use the mm10 genome.\n")
    cat("I am downloading the required annotation and sequence files. Please be patient.\n")
    download.file("http://steinmetzlab.embl.de/genomes/Mus_musculus.GRCm38.92.chr.gtf.gz", destfile = file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.92.chr.gtf.gz"))
    download.file("http://steinmetzlab.embl.de/genomes/Mus_musculus.GRCm38.73.dna.chromosome.all.spikes.fa.gz", destfile = file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.73.dna.chromosome.all.spikes.fa.gz"))
    download.file("http://steinmetzlab.embl.de/genomes/Mus_musculus.GRCm38.cdna.convenient.all.fa.gz", destfile = file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.cdna.convenient.all.fa.gz"))
    cat("Finished downloading, unpacking. Please be patient.\n")
    gunzip(file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.73.dna.chromosome.all.spikes.fa.gz"))
    gunzip(file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.cdna.convenient.all.fa.gz"))
    cat("Done downloading and unpacking.\n")
  }
  
  TxDb <- import.gff(file.path(script.basename,"genomeFiles", "Mus_musculus.GRCm38.92.chr.gtf.gz"),format="gtf")
  TxDb <- subset(TxDb, transcript_biotype == "protein_coding" & type == "exon") 
  
  org <- org.Mm.eg.db
  
  genome_fasta <- file.path(script.basename,"genomeFiles","Mus_musculus.GRCm38.73.dna.chromosome.all.spikes.fa")
  tx_fasta <- file.path(script.basename,"genomeFiles","Mus_musculus.GRCm38.cdna.convenient.all.fa")
  
} else if (genome == "hg19") {
  suppressMessages({
  pkgTest("BSgenome.Hsapiens.UCSC.hg19", bc =T)
  pkgTest("org.Hs.eg.db", bc =T)
  })
  Bsgenome <- BSgenome.Hsapiens.UCSC.hg19
  
  if(!file.exists(file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh37.87.chr.gtf.gz"))) {
    cat("This appears to be the first time that you use the hg19 genome.\n")
    cat("I am downloading the required annotation and sequence files. Please be patient.\n")
    download.file("http://steinmetzlab.embl.de/genomes/Homo_sapiens.GRCh37.87.chr.gtf.gz", destfile = file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh37.87.chr.gtf.gz"))
    download.file("http://steinmetzlab.embl.de/genomes/Homo_sapiens.GRCh37.68.dna.chromosomes.withERCC.fa.gz", destfile = file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh37.68.dna.chromosomes.withERCC.fa.gz"))
    download.file("http://steinmetzlab.embl.de/genomes/Human_ERCC_combined.convenient.fa.gz", destfile = file.path(script.basename,"genomeFiles", "Human_ERCC_combined.convenient.fa.gz"))
    cat("Finished downloading, unpacking. Please be patient.\n")
    gunzip(file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh37.68.dna.chromosomes.withERCC.fa.gz"))
    gunzip(file.path(script.basename,"genomeFiles", "Human_ERCC_combined.convenient.fa.gz"))
    cat("Done downloading and unpacking.\n")
  }
  
  TxDb <- import.gff(file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh37.87.chr.gtf.gz"),format="gtf")
  TxDb <- subset(TxDb, transcript_biotype == "protein_coding" & type == "exon") 
  org <- org.Hs.eg.db
  
  genome_fasta <- file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh37.68.dna.chromosomes.withERCC.fa")
  tx_fasta <- file.path(script.basename,"genomeFiles", "Human_ERCC_combined.convenient.fa")
} else  if (genome == "hg38") {
  suppressMessages({
    pkgTest("BSgenome.Hsapiens.UCSC.hg38", bc=T)
    pkgTest("org.Hs.eg.db", bc=T)
  })
  
  if(!file.exists(file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh38.89.chr.protein-coding.exon.gtf.gz"))) {
    cat("This appears to be the first time that you use the hg38 genome.\n")
    cat("I am downloading the required annotation and sequence files. Please be patient.\n")
    download.file("http://steinmetzlab.embl.de/genomes/Homo_sapiens.GRCh38.89.chr.protein-coding.exon.gtf.gz", destfile = file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh38.89.chr.protein-coding.exon.gtf.gz"))
    download.file("http://steinmetzlab.embl.de/genomes/hg38_ucsc.fasta.gz", destfile = file.path(script.basename,"genomeFiles", "hg38_ucsc.fasta.gz"))
    download.file("http://steinmetzlab.embl.de/genomes/Human_ERCC_combined.convenient.fa.gz", destfile = file.path(script.basename,"genomeFiles", "Human_ERCC_combined.convenient.fa.gz"))
    cat("Finished downloading, unpacking. Please be patient.\n")
    gunzip(file.path(script.basename,"genomeFiles", "hg38_ucsc.fasta.gz"))
    gunzip(file.path(script.basename,"genomeFiles", "Human_ERCC_combined.convenient.fa.gz"))
    cat("Done downloading and unpacking.\n")
  }
  
  TxDb <- import.gff(file.path(script.basename,"genomeFiles", "Homo_sapiens.GRCh38.89.chr.protein-coding.exon.gtf.gz"),format="gtf")
  TxDb <- subset(TxDb, transcript_biotype == "protein_coding" & type == "exon") 
  Bsgenome <- BSgenome.Hsapiens.UCSC.hg38
  org <- org.Hs.eg.db
  
  genome_fasta <- file.path(script.basename,"genomeFiles", "hg38_ucsc.fasta")
  tx_fasta <-file.path(script.basename,"genomeFiles", "Human_ERCC_combined.convenient.fa")
  
} else stop("Genome needs to be one of hg19, hg38 or mm10")

seqlevels(TxDb) <- paste0("chr",seqlevels(TxDb))
TxDb <- makeTxDbFromGRanges(TxDb)


optim <- !is.null(opt$optim) | !is.null(opt$optimThr)| !is.null(opt$optimThrSame)
optim.thr <- ifelse(!is.null(opt$optimThr), opt$optimThr, 15)
optim.thr.same <- ifelse(!is.null(opt$optimThrSame), opt$optimThrSame, 30)
optim.outer <- !is.null(opt$optimOuter)

blast <- !is.null(opt$blast) | !is.null(opt$txdb)

if (input == "cDNA") {
  if (!is.null(opt$txdb)) {
    blastdb_tx <- opt$txdb
  } else if (blast) {
    cat("Preparing transcriptome blast DB...\n")
    blastdb_tx <- make_blast_db(tx_fasta)
    cat("Created transcriptome blastdb for",genome,":",blastdb_tx,"\n")
  }
}


if (!is.null(opt$genomedb)) {
  blastdb_genome <- opt$genomedb
} else {
  blastdb_genome <- make_blast_db(genome_fasta)
  cat("Preparing genome blast DB...\n")
  cat("Created genome blastdb for",genome,":",blastdb_genome,"\n")
}


if ( is.null(opt$numprimers ) ) { numprimers <-  5 } else numprimers <- opt$numprimers

params_inner <- list()
params_outer <- list()
if (!is.null(opt$TmOuter)) params_outer$Tm <- as.numeric(strsplit(opt$TmOuter,",")[[1]]) else if (nested %in% c("RT","RTonly")) params_outer$Tm <- c(58,52,63) else params_outer$Tm <- c(60,55,65)
if (!is.null(opt$TmInner)) params_inner$Tm <- as.numeric(strsplit(opt$TmInner,",")[[1]]) else params_inner$Tm <- c(60,55,65)
if (!is.null(opt$nCharOuter)) params_outer$nChar <- as.integer(strsplit(opt$nCharOuter,",")[[1]]) else if (nested %in% c("RT","RTonly")) params_outer$nChar <- c(20,16,25) else params_outer$nChar <- c(23,19,27)
if (!is.null(opt$nCharInner)) params_inner$nChar <- as.integer(strsplit(opt$nCharInner,",")[[1]]) else params_inner$nChar <- c(20,16,25)



                          #in case RTonly, margin*mean(length_outer) is the minimal RT product length
#the script first designs the inner set and then it design the outer set around it.




###### 2. Prepare targets ######

#if the input is coordinates: if any target is within 50 bp of another target, make them a range
# convert to GRanges object
if (targetType == "coordinate") {
  targets <- targets2range(targets)
  
  #if the input is cDNA: extract all refseq transcripts containing the site
  if (input == "cDNA") {
    targets <- getcDNASeqfromCoordinate(targets, TxDb, Bsgenome,mean(length_inner),margin)
  } else {
     #otherwise, simply extract the sequencde surrounding the site
    targets <- getGenomicfromCoordinate(targets, Bsgenome,mean(length_inner),mean(length_outer),margin)
   }
  
} else {
  #if the input is a list of gene names: return targets object (i.e. Granges), let the position be a splice site.
  symbols <- targets
  targets <- getcDNASeqfromSymbol(symbols, mean(length_outer),margin, verbose = verbose) 
  
}
names(targets) <- sapply(targets, "[[", "id")

##### 3. Design inner primers ########
#now targets contains the sequence that we need to target together with a position.
#first, design inner primers surrounding that position.
if (nested != "RTonly") {
  
  cat("Running primer3 for cadidate identification (inner PCR)...\n")
  
  lines <- unlist(lapply(targets, getLines, numprimers = numprimers, range = length_inner, primerParams = params_inner))
  writeLines(lines, "pass2primer3.io")
  primer3 <- pipe(sprintf("%s pass2primer3.io | perl %s", file.path(primer3_path, "primer3_core"), file.path(script.basename,"parsePrimers.pl")))
  primer3result <- readLines(primer3)
  close(primer3)
  
  primer3result <- strsplit(primer3result, ";")
  for ( i in 1:length(primer3result)) {
    targets[[i]]$left <- strsplit(primer3result[[i]][2], ",")[[1]]
    targets[[i]]$right <- strsplit(primer3result[[i]][3], ",")[[1]]
  }
  
  cat("Done...\n")
  
  if (blast){
  cat("Blasting primers...\n")
  targets <- evaluate_blast_primers(targets, input, ifelse(input=="cDNA",blastdb_tx,blastdb_genome),verbose = verbose)
  } else {
    targets <- lapply(targets, function(x) {
      x$blast_inner_warn <- F
      x
    })
  }
  

  #identify for each primer pair if the left or right primer is with x bases of the mutation
  #add adapters for optimization run
  if (!is.null(opt$readlength)) {
    cat("Assuring that the first read / left primer is no more than", opt$readlength, "bases from the mutation\n")
    targets <- lapply(targets,addprefixes,firstreadprefix = prefixInnerLeft, secondreadprefix=prefixInnerRight, readlength = opt$readlength)
  } else {
    targets <- lapply(targets, function(x) {x$isleft <- rep(T,length(x$left)); x})
  }
  
  # save(targets, input, params_inner, file = "state.rda")
  # stop()
  if (optim) {
    cat("Finding ideal primer pairs (inner PCR):\n")
    cat("Checking", sum(sapply(targets,function(x) length(x$left)))*sum(sapply(targets,function(x) length(x$right))) , "possible combinations of primers for potential to form dimers, this may take a while. To speed up, decrease numprimers or switch off optimization")
    targets <- optimTargets(targets, input, primerParams = params_inner,verbose=verbose,firstreadprefix = prefixInnerLeft, secondreadprefix=prefixInnerRight, anyThr=optim.thr, anyThrSame = optim.thr.same)  
    cat("Done...\n")
  } else {
    targets <- lapply(targets, function(x) {
      x$left <- x$left[1]
      x$right <- x$right[1]
      x
    })
  }
  
  save(targets,file="inner.result.rda")
  
  
}
##### 4. Design outer primers ########


if (nested != "none" & nested != "RTonly") {
  #determine where in the sequence the selected inner primers bind.  
  targets <- lapply(targets, function(x) {
    us <- ifelse(nested == "RT", x$seq,x$completeSeq)
    a <- str_locate(us, gsub(prefixInnerLeft,"",x$left))[1,"start"]
    b <- str_locate(us, gsub(prefixInnerRight,"",x$right))[1,"start"]
    c <- str_locate(us,as.character(reverseComplement(DNAString(gsub(prefixInnerRight,"",x$right)))))[1,"end"]
    d <- str_locate(us,as.character(reverseComplement(DNAString(gsub(prefixInnerLeft,"",x$left)))))[1,"end"]
    x$end_inner_product <- max(c(a,b,c,d),na.rm = T)
    x$start_inner_product <- min(c(a,b,c,d),na.rm = T)
    
    if (x$end_inner_product > nchar(us) - 50) {
      c <- str_locate(us,as.character(reverseComplement(DNAString(gsub(prefixInnerRight,"",x$right)))))[1,"start"]
      d <- str_locate(us,as.character(reverseComplement(DNAString(gsub(prefixInnerLeft,"",x$left)))))[1,"start"]
      x$end_inner_product <- ifelse(is.na(c),d,c)
    }
    x
  })
}

if (nested == "RTonly") targets <- lapply(targets,function(x) {
  x$end_inner_product <- x$pos+1
  x$start_inner_product <- x$pos-1
  x$left <- NA
  x$right <- NA
  x$blast_inner_warn <- F
  x
  })
  
  #in case the outer primer is RT, it simply needs to be placed (anywhere) downstream of $pos
if (nested == "RT" | nested == "RTonly") {
  cat("Designing RT primers\n")
    lines <- unlist(lapply(targets, getLinesRT, numprimers = numprimers, range = length_outer, primerParams = params_outer))
    
    writeLines(lines, "pass2primer3.io")
    primer3 <- pipe(sprintf("%s pass2primer3.io | perl %s", file.path(primer3_path, "primer3_core"), file.path(script.basename,"parsePrimers.pl")))
    primer3result <- readLines(primer3)
    close(primer3)
    
    
    primer3result <- strsplit(primer3result, ";")
    for ( i in 1:length(primer3result)) {
      targets[[i]]$right_outer <-  strsplit(primer3result[[i]][ ifelse(targets[[i]]$sense,3,2 )], ",")[[1]]
      targets[[i]]$left_outer <- NA
    }
    
    #now blast, this time
    if (blast){
    cat("Blasting primers...\n")
      cat("Current default is to blast RT primers against the GENOME to avoid targeting non-coding transcripts...\n")
    targets <- evaluate_blast_primers_RT(targets, blastdb_genome,target="genome",verbose=verbose)
    } else {
      targets <- lapply(targets, function(x) {
        x$blast_outer_warn <- F
        x
      })
    } 
    
    if (optim) {
      cat("Finding ideal primers to avoid mispriming from other primers...:\n")
      cat("Checking", sum(sapply(targets,function(x) length(x$right_outer)))*sum(sapply(targets,function(x) length(x$right_outer))) , "possible combinations of primers for potential to form dimers, this may take a while. To speed up, decrease numprimers or switch off optimization")
      
      targets <- optimTargetsRT(targets, input,params_outer,verbose)
      } else {
        targets <- lapply(targets, function(x) {
          x$right_outer <- x$right_outer[1]
          x
          })
      }
  }

if (nested == "full") {
  cat("Running primer3 for cadidate identification (outer PCR)...\n")
  
  lines <- unlist(lapply(targets, getLinesOuter, numprimers = numprimers, range = length_outer, primerParams = params_outer))
  writeLines(lines, "pass2primer3.io")
  primer3 <- pipe(sprintf("%s pass2primer3.io | perl %s", file.path(primer3_path, "primer3_core"), file.path(script.basename,"parsePrimers.pl")))
  primer3result <- readLines(primer3)
  close(primer3)
  
  primer3result <- strsplit(primer3result, ";")
  for ( i in 1:length(primer3result)) {
    targets[[i]]$left_outer <- strsplit(primer3result[[i]][2], ",")[[1]]
    targets[[i]]$right_outer <- strsplit(primer3result[[i]][3], ",")[[1]]
  }
  
  cat("Done...\n")
  
  if (blast){
    cat("Blasting primers...\n")
    targets <- evaluate_blast_primers(targets, input, ifelse(input=="cDNA",blastdb_tx,blastdb_genome),mode="outer",verbose = verbose)
  } else {
    targets <- lapply(targets, function(x) {
      x$blast_outer_warn <- F
      x
    })
  }
  
  if (optim.outer) {
    cat("Finding ideal primer pairs (outer PCR):\n")
    cat("Checking", sum(sapply(targets,function(x) length(x$left_outer)))*sum(sapply(targets,function(x) length(x$right_outer))) , "possible combinations of primers for potential to form dimers, this may take a while. To speed up, decrease numprimers or switch off optimization")
    targets <- optimTargets(targets, input, primerParams = params_outer,mode="outer",verbose=verbose,anyThr = 30, anyThrSame = 30)  #outer primers: not as critical, use more relaxed thresholds.
    cat("Done...\n")
  } else {
    targets <- lapply(targets, function(x) {
      x$left_outer <- x$left_outer[1]
      x$right_outer <- x$right_outer[1]
      x
    })
  }
}

if (nested == "none") {
  targets <- lapply(targets, function(x) {
    x$right_outer <- NA
    x$left_outer <- NA
    x$blast_outer_warn <- F
    x
  })
}


output <- data.frame(
  id = sapply(targets, "[[","id"),
  left_inner = sapply(targets, "[[","left"),
  right_inner = sapply(targets, "[[","right"),
  left_outer = sapply(targets, "[[","left_outer"),
  right_outer = sapply(targets, "[[","right_outer"),
  blast_warning_inner =  sapply(targets, "[[","blast_inner_warn"),
  blast_warning_outer =  sapply(targets, "[[","blast_outer_warn"),
  exon_selection_warning = sapply(targets, "[[", "ExonSelectionWarning"),
  stringsAsFactors = F
)





#for output: blast all primers again, to define their genomic position and create a bed file.

if (!is.null(opt$bed)) {
  cat("Blasting final set of primners\n")
  # save.image("state3.rda")
  # stop("Saved")
  header <- "track name=\"PrimerDesign\" description=\"Primer Design\" visibility=2 itemRgb=\"On\""
  file.remove(opt$bed)
  writeLines(header,con=opt$bed)
  
  if(!all(is.na(output$left_inner))) {
    bed_left_inner <- final_blast_primers(output, "left_inner",blastdb_genome, prefixInnerLeft)
    write.table(bed_left_inner, file=opt$bed, quote=F,sep="\t",col.names = F, row.names=F,append=T)
    
  }
  if(!all(is.na(output$right_inner))) {
    bed_right_inner <- final_blast_primers(output, "right_inner",blastdb_genome, prefixInnerRight)
    write.table(bed_right_inner, file=opt$bed, quote=F,sep="\t",col.names = F, row.names=F,append=T)
    
  }
  if(!all(is.na(output$right_outer))) {
    bed_right_outer <- final_blast_primers(output, "right_outer",blastdb_genome)
    write.table(bed_right_outer, file=opt$bed, quote=F,sep="\t",col.names = F, row.names=F,append=T)
    
  }
  if(!all(is.na(output$left_outer))) {
    bed_left_outer <- final_blast_primers(output, "left_outer",blastdb_genome)
    write.table(bed_left_outer, file=opt$bed, quote=F,sep="\t",col.names = F, row.names=F,append=T)
    
  }
  
}

# 
# if(!is.null(opt$prefixInner)) {
#   if(!all(is.na(output$left_inner)))  output$left_inner <- paste0(opt$prefixInner,output$left_inner)
#   if(!all(is.na(output$right_inner))) output$right_inner <- paste0(opt$prefixInner,output$right_inner)
# }
# 
# if(!is.null(opt$prefixOuter)) {
#   if(!all(is.na(output$left_outer))) output$left_outer <- paste0(opt$prefixOuter,output$left_outer)
#   if(!all(is.na(output$right_outer))) output$right_outer <- paste0(opt$prefixOuter,output$right_outer)
# }

if (!is.null(opt$csv)) {
  write.csv(output,file = opt$csv,quote=F)
}
