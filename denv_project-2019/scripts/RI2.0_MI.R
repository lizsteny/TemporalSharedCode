### Libraries
library(methods)
library(plyr)
library(pbapply)
library(clue)
library(parallel)

options(bitmapType = 'cairo') # Make sure it works on ADA without support for X11

### Create list file with input data
files<- list(c("denv4/data/seqs/denv4_aa.fasta", # Fasta file
               0, # Reference position of sequences (leave 0 if not relevant)
               "TESTrun", # Output tag
               TRUE, TRUE, TRUE, TRUE, # Which subscripts to run (keep all true) 
               0.40 )) # Entropy cutoff before calculating MI 


### Create output folders
allelesOutFolder <- "denv4/data/output_alleles/"
MIOutFolder <- "denv4/data/output_MI/"
dir.create(file.path("denv4/data/output_alleles/"), showWarnings = FALSE)
dir.create(file.path("denv4/data/output_MI/"), showWarnings = FALSE)

### Define threholds
threholdX <- 34 # Convert loci <= cutoff to X!X
thresholdXs <- 0.4 # Reject when X!X >= proportion of samples (site of pairs)
thresholdAlleles <- 0.00 # Reject alleles below this proportion per locus
lociSeparator <- '!' # Separates loci alleles in output and calculations
alleleCollapse <- '+' # When making allele matrices squared, collapse alleles with this
gapToChar <- 'X' # All gaps '-' will be replaced by 'X' for convenience in interpretation
unwantedStates <- c('*','?') # There will be immediately converted to gap '-'

######

for(Lffd in files){

  ffd <- Lffd[1] # File to open
  reference <- Lffd[2] # Reference start position of sequences
  tag <- Lffd[3] # Tag to go on output files
  runGetAlleles <- Lffd[4] # Run script that gets the alleles
  runDigestAlleles <- Lffd[5] # Run the script that digests alleles
  runCalcIndexes <- Lffd[6] # Run the script that calculates MI for alleles
  runPlotNewMIPSallele <- Lffd[7] # Plot MI
  entropyTrheshold <- as.numeric(Lffd[8]) # Entropy cutoff before calculating MI
  
  ## Convert fasta to matrix format need for this script
  print(paste("<opening and converting FASTA file into a char matrix>", ffd))
  
  ## Giving a fake # separator allows to have data with line1=id, line2=seq, line3=id ...
  data <- read.table(ffd, 
                     header = FALSE, 
                     sep = "#", 
                     stringsAsFactors = FALSE, 
                     colClasses = "character", 
                     na.strings = "")
  
  ## Get sample names
  nameIndMatches <- apply(data, MARG = 1, FUN = function(x){grepl('>',x)})
  dataNames <- data[which(nameIndMatches),]
  
  ## Cicle all data lines and concatenate the sequence bits of each sample, which should be in between rows that we know are sample names
    prev <- TRUE
    allSeqs <- c()
    seq <- c()
    for(ii in 1:length(nameIndMatches)){
      curr<- nameIndMatches[ii]
      # if(curr==TRUE & prev==TRUE) -> do nothing (first step)
      if(curr == FALSE) seq <- paste0(seq, data[ii,], collapse = '', sep = '')
      if(curr == TRUE & prev == FALSE) {allSeqs <- rbind(allSeqs, seq); seq <- c()}
      if(ii == length(nameIndMatches)) allSeqs <- rbind(allSeqs, seq)
      prev <- curr
    }
    ## Separate / explode strings into each site / position
    X <- lapply(allSeqs, FUN = function(x){substring(x, seq(1,nchar(x),1), seq(1,nchar(x),1))})
    data <- matrix(unlist(X), nrow = length(X), byrow = TRUE)
    originaldataNaming <- (1:ncol(data))+as.numeric(reference) # Code will ignore the original names

######
    
  print("<converting unwanted states to gaps>")

   for(uw in unwantedStates){
     data[data == uw] <- '-'
   }

######

  miPlotThreshold <- 0.50 # Do not save output below this MI threshold
  savePropSchedule <- 0.1 # Save in temporary files every 100x X%
  reference <- as.numeric(reference)

  ##tag is built with all thresholds
  tag <- paste0(tag,'_Xs',thresholdXs,'_Al',thresholdAlleles,'_En',entropyTrheshold)

  if(runGetAlleles)  source("scripts/code_alleles/A1.6_allele_counting.R")

  if(runDigestAlleles) source("scripts/code_alleles/A2.4_alllele_counts_to_tables.R")

  if(runCalcIndexes) source("scripts/code_alleles/A3.6_calc_indexes_from_alleles.R")

  if(runPlotNewMIPSallele) source("scripts/code_alleles/A4.0_MIPS_allele.R")

}

## Clean GE
rm(list = ls(all.names = TRUE)) 