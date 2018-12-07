#check if installation is correct

if(!file.exists(file.path(script.basename,"config.txt")) {
  cat("No configuration file found.\n")
  cat("This may be this first time you run PrimerDesign.\n")
  cat("Please secify the path to the binaries folder of NCBI blast (e.g. /home/User/ncbi-blast-2.6.0/bin)\n")
  cat("BLAST can be downloaded from https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download\n")
  NCBI_path <- readLines(n=1)
  if(!file.exists(file.path(NCBI_path, "makeblastdb"))) stop("No file named makeblastdb in",NCBI_path,"\nPlease specify a correct path to NCBI blast binaries")
  
  cat("Please secify the path to the folder containing the primer3_core executable (e.g. /home/User/primer3/)\n")
  cat("Folder also needs to contain a subfolder named primer3_config; this should usually be the case.\n")
  cat("Primer3 and the config folder can be downloaded from https://github.com/primer3-org/primer3\n")
  primer3_path <- readLines(n=1)
  if(!file.exists(file.path(primer3_path, "primer3_core"))) stop("No file named primer3_core in", primer3_path,"\nPlease specify a correct path to primer3 binary")
  if(!file.exists(file.path(primer3_path,"primer3_config/"))) stop("Error - please make sure that the primer3 path contains a folder primer3_config. If for some reason this is not the case, download the folder from https://github.com/primer3-org/primer3/tree/master/src")
  
  writeLines(c(NCBI_path,primer3_path), con = file.path(script.basename,"config.txt"))
  
  if (tryCatch(system("perl -e \"print 'hallo'\"",intern = T), error="nope") =="nope" ) stop("Failed to run perl. Please make sure that the command perl is available in your environment.")
  
}) else {
  NCBI_path <- readLines(file.path(script.basename,"config.txt"))[1]
  primer3_path <- readLines(file.path(script.basename,"config.txt"))[2]
  if(!file.exists(file.path(NCBI_path, "makeblastdb")) | !file.exists(file.path(primer3_path, "primer3_core"))) {
    file.remove(file.path(script.basename,"config.txt"))
    stop("The config.txt file was incorrect. I have removed the file, please restart the program and specify correct paths to NCBI Blast and primer3.")
  }
  
}


pkgTest <- function(x, bc =F) {
  if (!require(x,character.only = TRUE))
  {
    cat("Package",x,"is missing, installing...")
    if (bc) biocLite(x) else install.packages(x,dep=TRUE)
    
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

source("https://bioconductor.org/biocLite.R")

pkgTest("getopt")
pkgTest("R.utils")
pkgTest("utils")
pkgTest("magrittr")
pkgTest("stringr")
pkgTest("GenomicRanges", bc=T) #bionductor
pkgTest("rtracklayer", bc =T) #bioconductor
pkgTest("GenomicFeatures", bc=T) #bioconductor




