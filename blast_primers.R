# by Andreas
## FUN: blast primer sequences against genome sequence and return found hits
#  for each primer.

#/g/steinmetz/gschwind/software/blast+/ncbi-blast-2.6.0+/bin/makeblastdb -in /g/steinmetz/genome/Homo_sapiens/GRCh38/fasta/GRCh38.primary_assembly.genome.ERCC.fa -input_type fasta -dbtype nucl -out /Volumes/steinmetz/velten/Scripts/Enhancer/blast-170630165622/blastdb
make_blast_db <- function(fasta, tmpdir = getwd(),
                          blast_root = "/g/steinmetz/gschwind/software/blast+/ncbi-blast-2.6.0+/bin") {
  tmpdir <- paste0(tmpdir, "/blast-", format(Sys.time(), "%y%m%d%H%M%S"))
  dir.create(tmpdir)
  
  # set parameters for blast data base generation
  cat("Creating blast db...")
  blastdb <- paste0(tmpdir ,"/blastdb")
  db_exec <- paste0(blast_root, "/makeblastdb")
  db_args <- paste("-in", fasta, "-input_type fasta -dbtype nucl -out",
                   blastdb)
  
  # create blast database
  cat(db_exec, db_args)
  system2(command = db_exec, args = db_args)
  tmpdir
}

blast_primers <- function(primers, tmpdir,
                          blast_root = "/g/steinmetz/gschwind/software/blast+/ncbi-blast-2.6.0+/bin"){
  
  # create unique tmpdir
  blastdb <- paste0(tmpdir ,"/blastdb")
  
  # create DNAStringSet object with primer sequences
  primer_seqs <- DNAStringSet(primers)
  
  # save to fasta file
  primers_fasta <- paste0(tmpdir, "/primer_seqs.fasta")
  writeXStringSet(primer_seqs, filepath = primers_fasta)
  
  # set blastn parameters
  blast_exec <- paste0(blast_root, "/blastn")
  blast_out <- paste0(tmpdir, "/blast_out.txt")
  blast_args <- paste("-query", primers_fasta, "-out", blast_out, "-db", blastdb, "-task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 40")
  
  # run blastn for short query sequences
  system2(command = blast_exec, args = blast_args)
  # read blast ouptut
  hits <- read.table(blast_out, stringsAsFactor = FALSE)

  # return found hits
  return(hits)
  
}


filter_blast <- function(blasttable, sequences) {
  blasttable <- split(blasttable, blasttable$V1)
  blasttable <- mapply(function(x,s) subset(x, V3 > 0.9 & V8 >= nchar(s)-1),blasttable, sequences, SIMPLIFY = F) #needs matches at 3' end of primer and needs at least 90% identity
  blasttable <- do.call(rbind,blasttable)
}

evaluate_blast_primers <- function(targets, input, blastdb, mode = "inner", verbose =F) {
  #blast each set of of primers
  warn <- ifelse(mode == "inner", "blast_inner_warn", "blast_outer_warn")
  left <- ifelse(mode == "inner", "left", "left_outer")
  right <- ifelse(mode == "inner", "right", "right_outer")
  for (i in 1:length(targets)) {
    blast_left <- blast_primers(targets[[i]][[left]], blastdb) %>% filter_blast(targets[[i]]$left)
    blast_right <- blast_primers(targets[[i]][[right]], blastdb) %>% filter_blast(targets[[i]]$right)
    blast_both <- merge(blast_left, blast_right, by = c("V1","V2"),suffixes = c(".left",".right"))
    #cDNA: Remnove entries on gene of interest
    blast_both$V2 <- gsub("\\..+","",x = blast_both$V2)
    if (input=="cDNA") {
      blast_both$geneID <- AnnotationDbi::select(org, keys=blast_both$V2, columns="SYMBOL", keytype="ENSEMBL")[,"SYMBOL"]
      blast_both <- subset(blast_both, geneID != unique(gsub("--.+","",targets[[i]]$id)))
      remove <- unique(as.integer(gsub("Query_","",blast_both$V1)))
      
    } else {
      blast_both <- subset(blast_both, abs(V10.left - V10.right) < 3000)
      remove <- table(blast_both$V1)
      remove <- names(remove)[remove > 1]
      targets[[i]]$blast <- subset(blast_both, !V1 %in% remove)
      remove <- as.integer(gsub("Query_","",remove))
      
    }

    
    if (all(1:length(targets[[i]]$left) %in% remove)) {
     if (verbose) cat("BLAST", targets[[i]]$id, ": No good primers found, keeping all\n") 
      targets[[i]][[warn]] <- T
    } else if (length(remove)==0) {
      if (verbose) cat("BLAST",targets[[i]]$id, ": All primers are good\n")
      targets[[i]][[warn]] <- F
    } else {
      if (verbose) cat("BLAST",targets[[i]]$id, ": removing primers",remove,"\n")
      targets[[i]][[left]] <- targets[[i]][[left]][-remove]
      targets[[i]][[right]] <- targets[[i]][[right]][-remove]
      targets[[i]][[warn]] <- F
    }
  }
  
  targets

}


evaluate_blast_primers_RT <- function(targets, blastdb,target = "tx", verbose =F, nForceKeep = 3) {
  #blast each set of of primers
  for (i in 1:length(targets)) {
    blast_right <- blast_primers(targets[[i]]$right_outer, blastdb) %>% filter_blast(targets[[i]]$right_outer)
    #cDNA: Remnove entries on gene of interest
    if (target == "tx") {
    blast_right$V2 <- gsub("\\..+","",x = blast_right$V2)
    
      blast_right$geneID <- AnnotationDbi::select(org, keys=blast_right$V2, columns="SYMBOL", keytype="ENSEMBL")[,"SYMBOL"]
      blast_right <- subset(blast_right, geneID != unique(gsub("--.+","",targets[[i]]$id)) )
      remove <- unique(as.integer(gsub("Query_","",blast_right$V1)))
      
      if (all(1:length(targets[[i]]$right_outer) %in% remove)) {
        if (verbose) cat("BLAST", targets[[i]]$id, ": No good primers found, keeping all\n") 
        targets[[i]]$blast_outer_warn <- T
      } else if (length(remove)==0) {
        if (verbose) cat("BLAST",targets[[i]]$id, ": All primers are good\n")
        targets[[i]]$blast_outer_warn <- F
      } else {
        if (verbose) cat("BLAST",targets[[i]]$id, ": removing primers",remove,"\n")
        targets[[i]]$right_outer <- targets[[i]]$right_outer[-remove]
        targets[[i]]$blast_outer_warn <- F
      }

    } else {
      remove <- table(blast_right$V1)
      if (sum(remove > 1) == length(targets[[i]]$right_outer)) {
        targets[[i]]$blast_outer_warn <- T
        use <- names(remove)[remove <= minn(remove,nForceKeep )]
        use <- as.integer(gsub("Query_","",use))
        targets[[i]]$right_outer <- targets[[i]]$right_outer[use]
        if (verbose) cat("BLAST", targets[[i]]$id, ": No good primers found, keeping best",nForceKeep,"\n") 
        
        
      } else  if (sum(remove > 1 ) > 0) {
        remove <- names(remove)[remove > 1]
        remove <- as.integer(gsub("Query_","",remove))
        targets[[i]]$blast_outer_warn <- F
        targets[[i]]$right_outer <- targets[[i]]$right_outer[-remove]
        if (verbose) cat("BLAST",targets[[i]]$id, ": removing primers",remove,"\n")
        
        
      } else {
        targets[[i]]$blast_outer_warn <- F
        if (verbose) cat("BLAST",targets[[i]]$id, ": All primers are good\n")
        
      }
    }
    
    #as there will be no further refinement in the RT case: immediately spit out first primer
  } 
    
   
  targets
}



final_blast_primers <- function(output,field, blastdb) {
  colors = c("left_inner" = "175,0,0",
             "right_inner" = "0,175,0",
             "right_outer" = "0,0,200",
             "left_outer" = "255,0,0")
  #blast each set of of primers
    blasted <- blast_primers(output[,field], blastdb)
    blasted <- split(blasted, blasted$V1)
    blasted <- lapply(blasted,function(x) subset(x, V4 > nchar(output[as.integer(gsub("Query_","",V1)),field]) - 2  ))
    blasted <- lapply(blasted,function(x) {
      x$id <- output$id[as.integer(gsub("Query_","",x$V1))]
      x
      })
    
    blasted <- do.call(rbind,blasted)
    bed <- data.frame(chr = paste0("chr",blasted$V2), start = apply(blasted[,c("V9","V10")],1,min), end = apply(blasted[,c("V9","V10")],1,max),
                      name = blasted$id, score = round(blasted$V3), strand = ifelse(blasted$V9 > blasted$V10, "-", "+"),
                      thickStart = min(blasted[,c("V9","V10")]), thickEnd = max(blasted[,c("V9","V10")]),
                      itemRgb = colors[field])
    bed
                    
}