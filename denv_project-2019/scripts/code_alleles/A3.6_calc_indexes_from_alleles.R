require(clue)

buildPairFromMatrix <- function(mat, pos){
  alleles<- c()
  counts<- c()
  for(ir in 1:nrow(mat)){
      alleles<- c(alleles, paste(rownames(mat)[ir], colnames(mat), sep = lociSeparator))
      counts<- c(counts, mat[ir,])
  }
  counts<- as.numeric(counts)
  O <- order(counts, decreasing = TRUE)
  alleles <- alleles[O]
  counts <- counts[O]
  p <- c(as.character(pos), as.vector((t(matrix(c(alleles, counts),ncol = 2)))))
  return(p)
}

addIntoMatrix<- function(Y, X, inds, perCol, onames){
  if(perCol){
      ynames <- colnames(Y)
      for(ii in 1:length(inds)){
          Y[,inds[ii]] <- Y[,inds[ii]] + X[,ii]
          ynames[inds[ii]] <- paste(ynames[inds[ii]], onames[ncol(Y)+ii], sep = alleleCollapse)
      }
      colnames(Y) <- ynames
  }else{
    ynames <- rownames(Y)
    for(ii in 1:length(inds)){
        Y[inds[ii],]<- Y[inds[ii],] + X[ii,]
        ynames[inds[ii]] <- paste(ynames[inds[ii]], onames[nrow(Y)+ii], sep = alleleCollapse)
    }
    rownames(Y) <- ynames
  }
  return(Y)
}

addAssymIntoSymMatrix <- function(M){
    whichMaxInArray <- function(X){ return(which(X == max(X))[1]) }
    minD <- min(ncol(M), nrow(M))
    if(ncol(M)>nrow(M)){ # Cut cols
        L <- M[,1:minD] # Left symmetric
        R <- M[,(minD + 1):ncol(M)]  # Right assymetric
        if(length(dim(R)) == 0) R <- matrix(R,ncol = length((minD + 1):ncol(M)))
        adds<- as.numeric(apply(R, 2, whichMaxInArray))
        O <- addIntoMatrix(Y = L, X = R, inds = adds, perCol = TRUE, onames = colnames(M))
    }else{ # Cur rows
        U <- M[1:minD,] # Upper symmetric
        B <- matrix(M[(minD + 1):nrow(M),], nrow = length((minD + 1):nrow(M))) # Bottom assymetric
        adds <- as.numeric(apply(B, 1, whichMaxInArray))
        O <- addIntoMatrix(Y = U, X = B, inds = adds, perCol = FALSE, onames = rownames(M))
    }
    return(O)
}

alleles2Matrix<- function(pair){ # Separate counts and alleles
    pairM <- pair[3:length(pair)]
    pairM <- pairM[which(!is.na(pairM))]
    pairMalleles <- pairM[seq(1, length(pairM), 2)]
    pairMcounts <- pairM[seq(2, length(pairM), 2)]
    ## Order the data by counts
    O<- order(pairMalleles)
    pairMalleles <- pairMalleles[O]
    pairMcounts <- as.numeric(pairMcounts[O])
    ## Get actual alleles per loci 1 and 2 (new code uses a separator for loci)
    sepsLoci<- as.numeric(gregexpr(pattern = lociSeparator, as.character(pairMalleles)))
    pairMallelesA <- (unique(substr(pairMalleles, 1, sepsLoci-1)))
    pairMallelesB <- (unique(substr(pairMalleles, sepsLoci + 1, nchar(as.character(pairMalleles)))))
    ## Work on this function now
    M <- matrix(rep(0, length(pairMallelesA)*length(pairMallelesB)), ncol = length(pairMallelesA), nrow = length(pairMallelesB), byrow = FALSE)
    M <- t(M) # Apply(M, 1, as.numeric)
    colnames(M) <- pairMallelesB
    rownames(M) <- pairMallelesA
    for(ii in 1:length(pairMalleles)){
        aA <- substr(pairMalleles[ii], 1, sepsLoci[ii] - 1)
        aB <- substr(pairMalleles[ii], sepsLoci[ii] + 1, nchar(as.character(pairMalleles))[ii])
        M[aA,aB] <- as.numeric(pairMcounts[ii])
    }
    return(M)
}

optimise_trace <- function(mat){

### For testing 
## mat <- matrix(c(5,20,7,0,2,0,1,48), ncol = 4, byrow = TRUE)
## colnames(mat) <- c('a','v','l','t')
## rownames(mat) <- c('q','k')
## mat <- t(mat)
###
    A <- mat
    Nrow <- nrow(mat); Ncol<- ncol(mat); trans = FALSE;
    if(Nrow!=Ncol){
          if(Ncol < Nrow){ mat<- t(mat); trans = TRUE; }
          ss <- as.numeric(solve_LSAP(mat, maximum = TRUE))
          mat <- data.frame(mat)
          A<- mat[,c(ss,which(!((1:ncol(mat)) %in% ss)))]
          if(trans){ A <- t(A); }
          # if(trans){ M<-t(M); }
    }else{
          ss <- as.numeric(solve_LSAP(mat, maximum=TRUE))
          A <- mat[,ss]
    }
    return(A)
}

print("<calc_indexes_from_alleles>")

filename<-paste0(allelesOutFolder,tag,".extraInfo.Rdata")
load(filename) # Load saveEvents

allPairRejectN <- 0
allAlleleRejectN <- 0
allPairsN <- 0
actualSave <- 0

if(saveEvents > 0){
    for(Nsave in 1:saveEvents){

      filename<- paste0(allelesOutFolder, tag, "_allele_table.N", Nsave, ".csv")
      print(paste0("<file in> ", filename))
      allelesMatrix <- read.table(filename, header = FALSE, sep = ',', stringsAsFactors = FALSE, colClasses = "character", na.strings = "")
      print(paste0("<file in> nrow:",nrow(allelesMatrix),' ncol:',ncol(allelesMatrix)))

      ###
      ## Do we need to reject some based on proportion of Xs?
      ## Get a list of T/F to whether the pair has content > thresholdXs in Xs
      print(paste("<rejecting> reject pairs above Xs threshold? ..."))

      filterXs <- function(X){
          ## All the positions of each possible allele pair
          cc = seq(3,length(X),2)
          ## Get the counts of all allele pairs that have an x at least
          alleleCheck <- function(X, A){
            allele = as.character(A[X])
            freq = as.numeric(A[X + 1])
            if(length(grep(gapToChar, allele, ignore.case = TRUE)) > 0 ) return(freq)
            return(0)
          }
          allelesWXs <- sum(as.numeric(unlist(lapply(cc, alleleCheck, A = X))))
          allelesAll <- sum(as.numeric(X[cc + 1]), na.rm = TRUE)
          allelesWXzFreq <- allelesWXs/allelesAll
          # Mark for rejection is number of pairs with x is above threhsold
          if(allelesWXzFreq > thresholdXs) return(TRUE)
          return(FALSE)
      }
      rejectList <- apply(allelesMatrix, 1, filterXs)
      Nreject <- sum(rejectList)
      ## For a final print on totals rejected at the tend
      allPairRejectN <- allPairRejectN + Nreject
      allPairsN <- allPairsN + nrow(allelesMatrix)
      ## Do actual rejection
      print(paste("<rejecting>", Nreject, " pairs above Xs threshold (", round(100*Nreject/nrow(allelesMatrix), 3), "% )"))
      if(Nreject > 0) allelesMatrix <- allelesMatrix[!rejectList,]
      if(nrow(allelesMatrix) == 0) next;

      ### Reject alleles, not site pairs, which are below a certain threshold
      print(paste("<rejecting> reject alleles below proportion threshold? ..."))

      filterAlleleProportions<- function(X){
        alleleInds <- seq(3, length(X), 2)
        countInds <- seq(4, length(X), 2)
        counts <- as.numeric(X[countInds])
        props <- counts/sum(counts, na.rm = TRUE)
        rejectProps <- which(props <= thresholdAlleles)
        return(alleleInds[rejectProps])
      }
      rejectList <- apply(allelesMatrix, 1, filterAlleleProportions)

      if(length(rejectList) > 0){
        ## Actual rejection
        Nreject <- 0
        for(rr in 1:nrow(allelesMatrix)){
          reject <- rejectList[[rr]]
          allAlleleRejectN <- allAlleleRejectN + length(reject)
          Nreject <- Nreject + length(reject)
          if(length(reject) > 0){
            allelesMatrix[rr,reject] <- NA # Allele entry
            allelesMatrix[rr,reject + 1] <- NA # Count entry
          }
        }
        print(paste("<rejecting>",Nreject," alleles below proportion threshold"))
      }


      cl_MIs <- rep(NA, nrow(allelesMatrix))
      cl_MInorms <- rep(NA, nrow(allelesMatrix))
      cl_morphAs <- rep(NA, nrow(allelesMatrix))
      cl_morphBs <- rep(NA, nrow(allelesMatrix))

      W_MIs <- rep(NA, nrow(allelesMatrix))
      W_MInorms <- rep(NA, nrow(allelesMatrix))
      W_morphAs <- rep(NA, nrow(allelesMatrix))
      W_morphBs <- rep(NA, nrow(allelesMatrix))

      Dprime <- rep(NA, nrow(allelesMatrix))
      R2 <- rep(NA, nrow(allelesMatrix))


      print("<calculating 'MIPS'>")
      pbar <- txtProgressBar(min = 0, max = nrow(allelesMatrix), style = 3)
      for(rr in 1:nrow(allelesMatrix)){

                ## Get relevant info
                pair <- allelesMatrix[rr,]
                M <- alleles2Matrix(pair)

                if(nrow(M) <= 1 | ncol(M) <= 1){
                    ## Reject this pair with only one allele in one loci
                    next;
                    # cl_MIs[rr]<- NA
                    # cl_MInorms[rr]<- NA
                    # cl_morphAs[rr]<- nrow(M)
                    # cl_morphBs[rr]<- ncol(M)
                    # W_MIs[rr]<- NA
                    # W_MInorms[rr]<- NA
                    # W_morphAs[rr]<- nrow(M)
                    # W_morphBs[rr]<- ncol(M)
                    # Dprime[rr]<- NA
                    # R2[rr]<- NA
                }else{
                    ## Calc pop gen for this pair, this will be used below to figure allele name lengths
                    colTrace <- TRUE
                    if(nrow(M)>ncol(M)) colTrace <- FALSE
                    ## For debug
                    doneWalker <- FALSE
                    ## If asym, then use Walker's tracing, if not, use pair directly
                    if(ncol(M)!=nrow(M)){
                        ## Step 1-3 of Walker's algo
                        M <- optimise_trace(M)
                        ## Steps 4-5 (add outside of sym matrix keeping diag max)
                        M <- addAssymIntoSymMatrix(M)
                        ## Step 6 need to return M to 'pair' format
                        pair <- buildPairFromMatrix(M, pair[1:2])
                        ## For debug
                        # doneWalker <- TRUE
                    }
                    ## With the sym matrix, revert back to a list of alleles to calc MIs
                    pairPositions <- pair[1:2]
                    pairValues <- pair[3:length(pair)]; pairValues= pairValues[which(pairValues!='')]
                    aaNames <- pairValues[seq(1, length(pairValues), 2)]
                    aaCounts <- pairValues[seq(2, length(pairValues), 2)]
                    ## Calc probs/freq of AA
                    sepsLoci <- as.numeric(gregexpr(pattern = lociSeparator, aaNames))
                    allAAs <- data.frame(AA1 = substr(aaNames, 1, sepsLoci - 1),
                                        AA2 = substr(aaNames, sepsLoci + 1, nchar(aaNames)),
                                        count = as.numeric(aaCounts), stringsAsFactors = FALSE)
                    ## By default, the MI code is not expecting zero entries, remove them here
                    w_allAAs <- allAAs[which(allAAs$count!=0),]
                    ## Calculating the morphs of each site
                    W_morphAs[rr] <- length(unique(w_allAAs$AA1))
                    W_morphBs[rr] <- length(unique(w_allAAs$AA2))
                    ## Calculating components of MI 
                    allProbs <- ddply(.data=w_allAAs, .variables=.(AA1 , AA2), .parallel = FALSE,  extraD = w_allAAs,
                      function(x, extraD) {
                        # cat(paste(x$AA1,x$AA2));cat('\n')
                        # cat('------\n')
                            pAB <- x$count/sum(extraD$count)
                            pB <- sum(extraD$count[which(extraD$AA2==x$AA2)])/sum(extraD$count)
                            pA <- sum(extraD$count[which(extraD$AA1==x$AA1)])/sum(extraD$count)
                            return(c(pAB,pA,pB)) } )
                    colnames(allProbs)[3:5] <- c("pAB","pA","pB")
                    wMI <- sum( allProbs[["pAB"]]*log(allProbs[["pAB"]]/(allProbs[["pA"]]*allProbs[["pB"]])) )
                    W_MIs[rr] <- wMI
                    W_MInorms[rr] <- wMI/log(min(c(W_morphAs[rr],W_morphBs[rr])))


                    ## Use original matrix -> calculates classic MI values
                    ## Get relevant info
                    pair <- allelesMatrix[rr,]
                    # M <- alleles2Matrix(pair)
                    pairPositions <- pair[1:2]
                    pairValues <- pair[3:length(pair)]; pairValues= pairValues[which(pairValues!='')]
                    aaNames <- pairValues[seq(1, length(pairValues), 2)]
                    aaCounts <- pairValues[seq(2, length(pairValues), 2)]
                    # Calc probs/freq of AA
                    sepsLoci <- as.numeric(gregexpr(pattern=lociSeparator,aaNames))
                    c_allAAs <- data.frame(AA1=substr(aaNames,1,sepsLoci-1),
                                          AA2=substr(aaNames,sepsLoci+1,nchar(aaNames)),
                                          count=as.numeric(aaCounts),stringsAsFactors=FALSE)
                    ## Calculating the morphs of each site
                    cl_morphAs[rr]<- length(unique(c_allAAs$AA1))
                    cl_morphBs[rr]<- length(unique(c_allAAs$AA2))
                    ## Calculating components of MI ################
                    # c_allAAs<- c_allAAs[1:14,]
                    allProbs<- ddply(.data=c_allAAs, .variables=.(AA1 , AA2), .parallel=FALSE,  extraD=c_allAAs,
                      function(x, extraD) {
                            # cat(paste(x$AA1,x$AA2));cat('\n')
                            # cat('------\n')
                            pAB <- x$count/sum(extraD$count)
                            pB <- sum(extraD$count[which(extraD$AA2==x$AA2)])/sum(extraD$count)
                            pA <- sum(extraD$count[which(extraD$AA1==x$AA1)])/sum(extraD$count)
                            return(c(pAB,pA,pB)) } )
                    colnames(allProbs)[3:5] <- c("pAB","pA","pB")
                    clMI <- sum(allProbs[["pAB"]]*log(allProbs[["pAB"]]/(allProbs[["pA"]]*allProbs[["pB"]])))
                    cl_MIs[rr] <- clMI
                    cl_MInorms[rr] <- clMI/log(min(c(cl_morphAs[rr],cl_morphBs[rr])))


                    ## LD calculations for dimorphic pairs, note that we use here the data that was organized for the classic MI calculations

                    if(cl_morphAs[rr] == 2 & cl_morphBs[rr] == 2){
                      # pA <- allProbs[["pA"]]
                      # pB <- allProbs[["pB"]]
                      # pAB <- allProbs[["pAB"]]
                      # pa <- 1-pA
                      # pb <- 1-pB
                      # DAB<- pAB - pA * pB ##LD for all 4 allele combinations (indepedndently)
                      # a<- apply(cbind(pAB,pA,pB), MARG=1, FUN=function(X){X[1] - X[2] * X[3]})
                      # ## or LD for the loci, D = freq(AB).freq(ab) – freq(Ab).freq(aB).
                      # a<- pAB[1]*pAB[4]- pAB[2]*pAB[3]
                        
                      ## Calculating D for a particular combination of AA, the first one this is a problem if not looking at the other 3 possible combinations as it depends on frequencies, but not a problem if eventually you end up using R2, which corrects for that (see above for code that actual gets D for every of the 4 combinations)
                      pA <- allProbs[["pA"]][1]
                      pB <- allProbs[["pB"]][1]
                      pAB <- allProbs[["pAB"]][1]
                      pa <- 1-pA
                      pb <- 1-pB
                      DAB <- pAB - pA * pB
                      if(DAB>0){
                        Dmax<- min( (allProbs[["pA"]]*(1-allProbs[["pB"]]))[1], ((1-allProbs[["pA"]])*allProbs[["pB"]])[1] )
                      }else{
                        Dmax<- max( -((allProbs[["pA"]]*allProbs[["pB"]])[1]), - (((1-allProbs[["pA"]])*(1-allProbs[["pB"]]))[1]) )
                      }
                      Dprime[rr] <- DAB/Dmax
                      R2[rr] <- (DAB^2)/((pA*pa)*(pB*pb))
                    }else{
                      Dprime[rr] <- NA
                      R2[rr] <- NA
                    }
                }
                setTxtProgressBar(pbar, rr)
                ## For debug
                if(doneWalker){
                  print(paste(clMI,wMI))
                  line <- readline()
                }

      } ## Close FOR each pair
      close(pbar)

      output<- data.frame(
                          siteA = allelesMatrix[,1],
                          siteB = allelesMatrix[,2],
                          MI = cl_MIs,
                          MInorm = cl_MInorms,
                          morphA = cl_morphAs,
                          morphB = cl_morphBs,
                          wMI = W_MIs,
                          wMInorm = W_MInorms,
                          wmorphA = W_morphAs,
                          wmorphB = W_morphBs,
                          Dprime = Dprime,
                          R2 = R2
                          )

      ## At this point, many rows of output may be just NA, a consequence of having allocated vectors in the beginning with possible max size so, here, chop off rows that have not been filled in.
      subOutput <- output[,3:ncol(output)]
      isUnfilled <- function(x){
        return(sum(is.na(x)) == length(x))
      }
      accept <- !(apply(subOutput, MARG = 1, FUN = isUnfilled))
      output <- output[accept,]

      ## Keep appending output into a single table
      filename <- paste0(MIOutFolder, tag, "MI.csv")
      appFlag <- TRUE; if(actualSave == 0) appFlag <- FALSE
      write.table(output, filename, row.names = FALSE, col.names=!appFlag, sep = ',',append = appFlag)
      actualSave <- actualSave + 1
    }
}

print(paste("<rejected> total of ", allPairRejectN, "pairs due to Xs thresholds"))
print(paste("<rejected> total of ", allAlleleRejectN, "alleles due to proportion threshold"))
filename <- paste0(MIOutFolder, tag, "MI_rejected.csv")
write.table(data.frame(totalPairs = allPairsN, rejectedPairs = allPairRejectN, allAlleles = allAlleleRejectN), filename, col.names = TRUE, row.names = FALSE, sep = ',')
