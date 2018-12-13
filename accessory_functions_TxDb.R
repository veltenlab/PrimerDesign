maxn <- function(v, n) {
  v <- v[order(v,decreasing=T)]
  if(length(v) >=n ) v[n] else v[length(v)]
}

minn <- function(v, n) {
  v <- v[order(v,decreasing=F)]
  if(length(v) >=n ) v[n] else v[length(v)]
}


targets2range <- function(targets, range=50) {
  #Creates a Grange object, combines mutations within range bases into a single target
  targets <- strsplit(targets, ":")
  targets <- data.frame(V1 = sapply(targets, "[", 1),
                        V2 = as.integer(sapply(targets, "[", 2)))
  targets <- split(targets, targets$V1)
  targets <- lapply(targets, function(x) {
    if (nrow(x) == 1) data.frame(chr = unique(x$V1), start = x$V2, end = x$V2) else {
      x <- x[order(x$V2),]
      start <- x$V2[1]
      out_start <- c(); out_end <- c()
      for ( i in 2:nrow(x)) {
        if (x$V2[i] > start + range) {
          out_start <- c(out_start, start)
          out_end <- c(out_end, x$V2[i-1])
          start <- x$V2[i]
        } 
      }
      out_start <- c(out_start, start); out_end <- c(out_end, x$V2[i])
      data.frame(chr = unique(x$V1), start = out_start, end = out_end)
    }
  })
  targets <- do.call(rbind, targets)
  targets <- GRanges(seqnames = targets$chr, ranges = IRanges(targets$start,targets$end),"*")
  targets
}

getGenomicfromCoordinate <- function(targets, Bsgenome, rangeInner, rangeOuter, margin = 2) {
  r <- rangeInner * margin
  targets$seq <-as.character(getSeq(Bsgenome, targets+r))
  targets$pos <- r
  targets$endpos <- r + width(targets)
  targets$id <- as.character(targets)
  
   r <- ifelse(is.na(rangeOuter), 5000, max(c(5000,rangeOuter * margin)))
   targets$completeSeq <- as.character(getSeq(Bsgenome, targets+r))
   out <- lapply(targets, function(x) {
     list(Range = reduce(x),
          seq = x$seq,
          pos = x$pos,
          endpos = x$endpos,
          completeSeq = x$completeSeq,
          ExonSelectionWarning = F,
          id = x$id)
   })
}

getcDNASeqfromCoordinate <- function(targets, TxDb, Bsgenome, rangeInner, margin = 2) {
  #1. identify exon with mutation
  #2. identify transcripts containing that exons
  #3. select all exons common to these transcripts
  #--> In case of exon skipping, the script possibly combines exons that aren't directly linked but I never observed this behavior
  cds <- exons(TxDb, columns = c("tx_id","exon_id","gene_id"))
  exonTargets <- findOverlaps(targets, cds)
  
  remove <- which(!(1:length(targets) %in% queryHits(exonTargets)))
  if (length(remove) > 0) {
    warning("The following targets are not on exons and will be removed:", paste(remove, collapse = ", "))
  targets <- targets[unique(queryHits(exonTargets))]
  exonTargets <- findOverlaps(targets, cds)
  }
  
  #1. for all transcripts, get all exons
  exonTargets <- split(exonTargets, queryHits(exonTargets))
  sequences <- lapply(exonTargets, function(x) {
    #identify all exons associated with the transcript
    if (length(unique(unlist(cds[subjectHits(x)]$gene_id))) > 1) {
      cat("WARNING: ",as.character(targets[unique(queryHits(x))]), "Overlaps multiple genes:\n", paste(unique(unlist(cds[subjectHits(x)]$gene_id)),sep=","), "\nPlease enter gene id of the one to select:\n")
      use <- readLines(file("stdin"),n=1)
      x <- x[unlist(cds[subjectHits(x)]$gene_id) == use]
    }
    
    sense <- unique(strand(cds[subjectHits(x)])) == "+"
    
    
    txids <- lapply(subjectHits(x), function(y) cds$tx_id[[y]])
    txids <- txids[[which.max(sapply(txids,length))]]
    thisexons <- exons(TxDb, columns = c("tx_id","exon_id","gene_id"), filter = list(tx_id = txids))
    if (!width(thisexons[sapply(thisexons$tx_id, function(y) all(txids %in% y))]) < rangeInner * margin) thisexons <- thisexons[sapply(thisexons$tx_id, function(y) all(txids %in% y))] #exons that occur in all transcripts only
    
    genes <- unique(unlist(thisexons$gene_id))
    symbol <- AnnotationDbi::select(org, keys=genes, columns="SYMBOL", keytype="ENSEMBL")[,"SYMBOL"]
    thisexons <- reduce(thisexons)
    target <- targets[unique(queryHits(x))]
    whoistargetExon <- subjectHits(findOverlaps(target, thisexons))
    if (length(whoistargetExon)>1) stop("Something wromng, more than 1 CDS hit?")
    targetExon <- thisexons[whoistargetExon]
    if (whoistargetExon==1) tostart <- 0 else tostart <- sum(width(thisexons[1:(whoistargetExon-1)])) 
    #if (whoistargetExon==length(thisexons)) fromend <- 0 else fromend <- sum(width(thisexons[(whoistargetExon+1):length(thisexons)])) 
    #get 300 bp surrounding hit
    fromstart <- start(target) - start(targetExon)
    #toend <- end(targetExon) - end(target)
    sequences <- getSeq(Bsgenome, thisexons)
    #all information up to know was in genomic coordinate, but getSeq returns the true sequence
    #of course sequences are ordered iby gernomic coordinate but each exon is reverse complement
    #relative to + strand if the gene is on the minus strand
    if (unique(strand(thisexons))== "-") sequences <- reverseComplement(sequences) #<<--- double check behavior!
    sequences <- paste(sequences, collapse="")
    r <- rangeInner * margin
    startpos <- tostart+fromstart-r+1
    endpos <- tostart+fromstart+width(target)+r
    innerSeq <- substr(sequences,startpos,endpos)
    posInString <- ifelse(startpos > 0, r, tostart+fromstart)
    list(innerSeq,posInString,sequences,symbol,sense)
  })
  
  targets$seq <- sapply(sequences, "[[",1)
  targets$pos <- sapply(sequences, "[[",2)
  targets$endpos <- sapply(sequences, "[[",2) #check how to go about range targets
  targets$id <- paste(sapply(sequences, "[[",4), as.character(targets),sep="--")
  targets$completeSeq <- sapply(sequences, "[[",3)
  targets$sense <- sapply(sequences, "[[",5)
  
  
  out <- lapply(targets, function(x) {
    list(Range = reduce(x),
         seq = x$seq,
         pos = x$pos,
         endpos = x$endpos,
         completeSeq = x$completeSeq,
         ExonSelectionWarning = F,
         id = x$id,
         sense = x$sense)
  })
}

getcDNASeqfromSymbol <- function(symbols, rangeOuter, margin = 2,verbose =F) {
  geneIDs <- AnnotationDbi::select(org, keys=symbols, columns="ENSEMBL", keytype="SYMBOL")
  cds <- exons(TxDb, columns = c("tx_id","exon_id","gene_id"), filter = list(gene_id = geneIDs$ENSEMBL))
  cds$symbol <- AnnotationDbi::select(org, keys=unlist(cds$gene_id), columns="SYMBOL", keytype="ENSEMBL")[,2]
  
  cds <- split(cds, cds$symbol) #
  targets <- lapply(cds, function(x) {
    warn <- F
    transcripts <- unique(unlist(x$tx_id))
    use <- sapply(x$tx_id, function(y) all(transcripts %in% y))
    if (sum(use == 0)) {
      use <- max(sapply(x$tx_id, length)) == sapply(x$tx_id, length)
      usexons <- c() 
    } else {
    #get all adjacent exons
    usexons <- union( which(c(F,use) & c(use,F)) - 1,
                      which(c(F,use) & c(use,F))
                         ) #is the next also in all transcripts?
    }
  
    if (unique(strand(x)) == "+") {
      #plus stranded transcript: get the last two exons adjacent to each other
      #if useexons is  emtpy: instead use the exon with the highest number of txs
      if (length(usexons) > 0) target <- x[usexons]  else target <- x[maxn(which(use),2)]
      #get sequence
      seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, target) %>% paste(collapse="")
      #completeSeq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, x) %>% paste(collapse="")
      #note splice site
      if (length(target)>1 ) splicesite <- nchar(seq) - width(target)[length(target)] else splicesite <- NA

    } else {
      #minus stranded transcript: get the first two exons adjacent to each other
      if (length(usexons) > 0) target <- x[usexons]  else target <- x[minn(which(use),2)]
      seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, target) %>% rev %>% paste(collapse="")
      #completeSeq <- reverseComplement(getSeq(BSgenome.Mmusculus.UCSC.mm10, x)) %>% paste(collapse="")
      #note splice site
      if (length(target) > 1) splicesite <- nchar(seq) - width(target)[1] else splicesite <- NA

    }
    lengths <- sapply(x$tx_id, length)
      if (!is.na(rangeOuter)) {
        if (nchar(seq) < margin * rangeOuter )  {
          if (verbose) cat("For", unique(unlist(x$symbol)), ": Ideal exons have insufficient lengths, using exons surrounding ideal exon\n")
          warn <- T
          #reduce, then take neighboring exons and set splice site smartly
          reduced <- reduce(x)
          if (unique(strand(x)) == "-") reduced <- rev(reduced)
          
          
            isTarget <- max(which(overlapsAny(reduced, target)))
              cumlength <-  cumsum(width(reduced))
              if (cumlength[isTarget] > rangeOuter * margin) {
                #is sequence upstream of target is of sufficient length, place RT primer in target
                target <- reduced[1:isTarget]
                splicesite <- cumlength[isTarget-1]
              } else {
                #not so ideal case - RT primer cannot be placed in target
                if (verbose)  cat("For", unique(unlist(x$symbol)), ":RT primer cannot be placed in ideal exon\n")
                
                if (max(cumlength) < rangeOuter * margin) {
                  target <- reduced 
                  splicesite <- cumlength[length(cumlength)-1]
                  } else {
                  target <- reduced[1:min(which(cumlength > rangeOuter * margin))]
                  splicesite <- cumlength[min(which(cumlength > rangeOuter * margin))-1]
                }
              }
              seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, target) %>% paste(collapse="")

          
        } else if (verbose) {
          cat("For", unique(unlist(x$symbol)), ": Ideal case\n")
          
        }
      }
      

    list(seq = seq, pos = splicesite, endpos = splicesite, id = unique(x$symbol), completeSeq = seq, ExonSelectionWarning = warn, sense = unique(strand(x)) == "+")
    #list(seq = seq, completeSeq=seq, splicesite= splicesite, gene = unique(x$symbol), exons = target)
  })
  
}



getLinesRT <- function(x, numprimers,range, primerParams) {
  c(sprintf("SEQUENCE_ID=%s",x$id),
    sprintf("SEQUENCE_TEMPLATE=%s",ifelse(x$sense, substr(x$seq,x$end_inner_product,nchar(x$seq)),
            substr(x$seq,1,x$start_inner_product))),
    "PRIMER_TASK=generic",
    sprintf("PRIMER_PICK_LEFT_PRIMER=%d",ifelse(x$sense, 0,1)),
    "PRIMER_PICK_INTERNAL_OLIGO=0",
    sprintf("PRIMER_PICK_RIGHT_PRIMER=%d",ifelse(x$sense, 1,0)),
    sprintf("PRIMER_OPT_SIZE=%d",primerParams$nChar[1]), #was 3 bases shorter
    sprintf("PRIMER_MIN_SIZE=%d",primerParams$nChar[2]),
    sprintf("PRIMER_MAX_SIZE=%d",primerParams$nChar[3]),
    sprintf("PRIMER_OPT_TM=%f",primerParams$Tm[1]),
    sprintf("PRIMER_MIN_TM=%f",primerParams$Tm[2]),
    sprintf("PRIMER_MAX_TM=%f",primerParams$Tm[3]),
    sprintf("PRIMER_PRODUCT_SIZE_RANGE=%d-%d",range[1],range[2]),
    "P3_FILE_FLAG=0",
    "PRIMER_EXPLAIN_FLAG=1",
    sprintf("PRIMER_NUM_RETURN=%d",numprimers),
    sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s",file.path(primer3_path,"primer3_config/")),
    "=")
}



getLines <- function(x, numprimers, range, primerParams) {
  c(sprintf("SEQUENCE_ID=%s",x$id),
    sprintf("SEQUENCE_TEMPLATE=%s",x$seq),
    sprintf("SEQUENCE_TARGET=%d,%d",x$pos, x$endpos-x$pos+3),
    sprintf("SEQUENCE_EXCLUDED_REGION=%d,%d",x$pos-3,  x$endpos-x$pos+6),
    
    "PRIMER_TASK=generic",
    "PRIMER_PICK_LEFT_PRIMER=1",
    "PRIMER_PICK_INTERNAL_OLIGO=0",
    "PRIMER_PICK_RIGHT_PRIMER=1",
    sprintf("PRIMER_OPT_SIZE=%d",primerParams$nChar[1]), #was 3 bases shorter
    sprintf("PRIMER_MIN_SIZE=%d",primerParams$nChar[2]),
    sprintf("PRIMER_MAX_SIZE=%d",primerParams$nChar[3]),
    sprintf("PRIMER_OPT_TM=%f",primerParams$Tm[1]),
    sprintf("PRIMER_MIN_TM=%f",primerParams$Tm[2]),
    sprintf("PRIMER_MAX_TM=%f",primerParams$Tm[3]),
    sprintf("PRIMER_PRODUCT_SIZE_RANGE=%d-%d", range[1],range[2]),
    "P3_FILE_FLAG=0",
    "PRIMER_EXPLAIN_FLAG=1",
    sprintf("PRIMER_NUM_RETURN=%d",numprimers),
    sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s",file.path(primer3_path,"primer3_config/")),
    "=")
}

getLinesOuter <- function(x, numprimers, range, primerParams) {
  c(sprintf("SEQUENCE_ID=%s",x$id),
    sprintf("SEQUENCE_TEMPLATE=%s",x$completeSeq),
    sprintf("SEQUENCE_TARGET=%d,%d",x$start_inner_product, x$end_inner_product-x$start_inner_product+3),
    sprintf("SEQUENCE_EXCLUDED_REGION=%d,%d",x$start_inner_product, x$end_inner_product-x$start_inner_product),
        "PRIMER_TASK=generic",
    "PRIMER_PICK_LEFT_PRIMER=1",
    "PRIMER_PICK_INTERNAL_OLIGO=0",
    "PRIMER_PICK_RIGHT_PRIMER=1",
    sprintf("PRIMER_OPT_SIZE=%d",primerParams$nChar[1]), #was 3 bases shorter
    sprintf("PRIMER_MIN_SIZE=%d",primerParams$nChar[2]),
    sprintf("PRIMER_MAX_SIZE=%d",primerParams$nChar[3]),
    sprintf("PRIMER_OPT_TM=%f",primerParams$Tm[1]),
    sprintf("PRIMER_MIN_TM=%f",primerParams$Tm[2]),
    sprintf("PRIMER_MAX_TM=%f",primerParams$Tm[3]),
    sprintf("PRIMER_PRODUCT_SIZE_RANGE=%d-%d", range[1],range[2]),
    "P3_FILE_FLAG=0",
    "PRIMER_EXPLAIN_FLAG=1",
    sprintf("PRIMER_NUM_RETURN=%d",numprimers),
    sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s",file.path(primer3_path,"primer3_config/")),
    "=")
}



getLinesPairs <- function(id_fwd,id_rev,seq_fwd,seq_rev, primerParams) {
  c(sprintf("SEQUENCE_ID=fwd_%s-rev_%s", id_fwd, id_rev),
    sprintf("SEQUENCE_PRIMER=%s", seq_fwd),
    sprintf("SEQUENCE_PRIMER_REVCOMP=%s", seq_rev),
    "PRIMER_TASK=check_primers",
    "PRIMER_EXPLAIN_FLAG=1",
    sprintf("PRIMER_OPT_SIZE=%d",primerParams$nChar[1]), #was 3 bases shorter
    sprintf("PRIMER_MIN_SIZE=%d",primerParams$nChar[2]),
    sprintf("PRIMER_MAX_SIZE=%d",primerParams$nChar[3]),
    sprintf("PRIMER_OPT_TM=%f",primerParams$Tm[1]),
    sprintf("PRIMER_MIN_TM=%f",primerParams$Tm[2]),
    sprintf("PRIMER_MAX_TM=%f",primerParams$Tm[3]),
    sprintf("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s",file.path(primer3_path,"primer3_config/")),
    "=")
}


optimTargets <- function(targets, input, primerParams,mode="inner",verbose=F) {
 left <- ifelse(mode == "inner","left","left_outer")
 right <- ifelse(mode == "inner","right","right_outer")
  #could add blasting step here
  
  targets <- targets[sapply(targets, function(x) length(x[[left]])) > 0]
  #create a named list of all left and right primers
  allleft <- unlist(lapply(targets, function(x) {
    out <- x[[left]];
    names(out) <- paste(x$id, 1:length(out), sep="_")
    out
  }))
  allright <-  unlist(lapply(targets, function(x) {
    out <- x[[right]];
    names(out) <- paste(x$id, 1:length(out), sep="_")
    out
  }))
  
  #Check that no left primer is the reverse complement of a right primer!
  require(Biostrings)
  right_revcomp <- sapply(allright, function(x) as.character(reverseComplement(DNAString(x))))
  remove <- which(right_revcomp %in% allleft)
  or_remove <- which(allleft %in% right_revcomp)
  if (length(or_remove) < length(remove)) remove <- or_remove
  
  if (length(remove) > 0) {
    allright <- allright[-remove]
    allleft <- allleft[-remove]
  }
  
  #lines <- list()
  file.remove("pass2primer3_chk.io")
  outwriter <- file("pass2primer3_chk.io", open="a")
  #writeLines(unlist(lines), "pass2primer3_chk.io")
  for ( i in 1:length(allleft)) {
    for (j in 1:length(allright)) {
      #lines[[length(lines)+1]]
      line <- getLinesPairs(names(allleft)[i], names(allright[j]),
                            allleft[i], allright[j], primerParams = primerParams)
      writeLines(line,con = outwriter)
    }
  }
  close(outwriter)
  
  cat("Running primer3 for dimer identification...\n")
  
  
  primer3 <- pipe(sprintf("%s pass2primer3_chk.io | perl %s", file.path(primer3_path, "primer3_core"), file.path(script.basename,"parsePrimerCheck.pl")))
  primer3result <- read.csv(primer3,header = F,sep=",")
  
  #set up a matrix of primer-to-primer interactions
  any <- matrix(primer3result$V2, nrow = length(allleft), ncol = length(allright), dimnames = list(names(allleft), names(allright)), byrow = T)
  end <- matrix(primer3result$V3, nrow = length(allleft), ncol = length(allright), dimnames = list(names(allleft), names(allright)), byrow = T)
  save(any,end,allleft, allright, targets,file="primer3result.rda")
  #try to identify a vector which selects a) 1 primer pair per target and b) minimizes the outputdim
  #iteritavely
  #throw out the primer pair with the worst performance
  
  protected <- rep(F, nrow(any))
  counter <- 1
  while (nrow(any) > length(targets)){
    cat("Eliminating primer... (round", counter,")\n")
    counter <- counter+1
    maxis_left <- apply(any[!protected,!protected],1,max)
    maxis_right <- apply(any[!protected,!protected],2,max)
    if (sum(maxis_left) == 0) break;
    which_maxi_left <- which.max(maxis_left)
    which_maxi_right <- which.max(maxis_right)
    
    total_left <- sum(any[!protected,!protected][which_maxi_left,])
    total_right <- sum(any[!protected,!protected][,which_maxi_right])
    
    
    
    if (total_left > total_right) {
      any <- any[rownames(any) !=  names(which_maxi_left) , colnames(any) !=  names(which_maxi_left)]
      if (verbose) cat("Discarding primer pair", names(which_maxi_left) , "\n")
    } else  {
      any <- any[rownames(any) !=  names(which_maxi_right) , colnames(any) !=  names(which_maxi_right)]
      if (verbose) cat("Discarding primer pair", names(which_maxi_right) , "\n")
    }
    
    #test if only a single pair is left for a given target
    t <- sapply(strsplit(rownames(any), ".", fixed = T),"[",1)
    arrived<-names(table(t))[table(t) == 1]
    protected <- t%in%arrived
    
  }
  
  #now, for each target, take the primer pairs with the minimal interaction score
  u <- c()
  for (target in unique(t)) {
    consider <- which(t == target) 
    if (length(consider) == 1) {
      if (verbose) print(rownames(any)[consider])
      u <- c(u,rownames(any)[consider])
      targets[[target]][[left]] <- allleft[rownames(any)[consider]]
      targets[[target]][[right]] <- allright[rownames(any)[consider]]
    } else {
      use <- which.min(apply(any[consider,],1,sum))
      if (verbose) print(rownames(any)[consider][use])
      u <- c(u,rownames(any)[consider][use])
      targets[[target]][[left]] <- allleft[rownames(any)[consider][use]]
      targets[[target]][[right]] <- allright[rownames(any)[consider][use]]
    }
  }
  targets
}



optimTargetsRT <- function(targets, input, primerParams,verbose=F) {
  #could add blasting step here
  
  targets <- targets[sapply(targets, function(x) length(x$right_outer)) > 0]
 
  allright <-  unlist(lapply(targets, function(x) {
    out <- x$right_outer;
    names(out) <- paste(x$id, 1:length(out), sep="_")
    out
  }))
 
  file.remove("pass2primer3_chk.io")
  outwriter <- file("pass2primer3_chk.io", open="a")
  #writeLines(unlist(lines), "pass2primer3_chk.io")
  for ( i in 1:length(allright)) {
    for (j in 1:length(allright)) {
      #lines[[length(lines)+1]]=
      line <- getLinesPairs(names(allright)[i], names(allright[j]),
                            allright[i], allright[j], primerParams = primerParams)
      writeLines(line,con = outwriter)

    }
  }
  close(outwriter)
  
  cat("Running primer3 for dimer identification...\n")
  
  primer3 <- pipe(sprintf("%s pass2primer3_chk.io | perl %s", file.path(primer3_path, "primer3_core"), file.path(script.basename,"parsePrimerCheck.pl")))
  primer3result <- read.csv(primer3,header = F,sep=",")
  
  #set up a matrix of primer-to-primer interactions
  any <- matrix(primer3result$V2, nrow = length(allright), ncol = length(allright), dimnames = list(names(allright), names(allright)), byrow = T)
  end <- matrix(primer3result$V3, nrow = length(allright), ncol = length(allright), dimnames = list(names(allright), names(allright)), byrow = T)
  save(any,end, allright, targets,file="allright.rda")
  #try to identify a vector which selects a) 1 primer pair per target and b) minimizes the outputdim
  #iteritavely
  #throw out the primer pair with the worst performance
  
  protected <- rep(F, nrow(any))
  counter <- 1
  while (nrow(any) > length(targets)){
    cat("Eliminating primer... (round", counter,")\n")
    counter <- counter+1
    maxis_left <- apply(any[!protected,!protected],1,max)
    maxis_right <- apply(any[!protected,!protected],2,max)
    if (sum(maxis_left) == 0) break;
    which_maxi_left <- which.max(maxis_left)
    which_maxi_right <- which.max(maxis_right)
    
    total_left <- sum(any[!protected,!protected][which_maxi_left,])
    total_right <- sum(any[!protected,!protected][,which_maxi_right])
    
    
    
    if (total_left > total_right) {
      any <- any[rownames(any) !=  names(which_maxi_left) , colnames(any) !=  names(which_maxi_left)]
      if (verbose) cat("Discarding primer pair", names(which_maxi_left) , "\n")
    } else  {
      any <- any[rownames(any) !=  names(which_maxi_right) , colnames(any) !=  names(which_maxi_right)]
      if (verbose) cat("Discarding primer pair", names(which_maxi_right) , "\n")
    }
    
    #test if only a single pair is left for a given target
    t <- sapply(strsplit(rownames(any), ".", fixed = T),"[",1)
    arrived<-names(table(t))[table(t) == 1]
    protected <- t%in%arrived
    
  }
  
  #now, for each target, take the primer pairs with the minimal interaction score
  u <- c()
  for (target in unique(t)) {
    consider <- which(t == target) 
    if (length(consider) == 1) {
      if (verbose) print(rownames(any)[consider])
      u <- c(u,rownames(any)[consider])
      targets[[target]]$right_outer <- allright[rownames(any)[consider]]
    } else {
      use <- which.min(apply(any[consider,],1,sum))
      if (verbose) print(rownames(any)[consider][use])
      u <- c(u,rownames(any)[consider][use])
      targets[[target]]$right_outer <- allright[rownames(any)[consider][use]]
    }
  }
  targets
}