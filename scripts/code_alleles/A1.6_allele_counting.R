## Download from: http://winvector.github.io/Parallel/bindToEnv.R
source("scripts/code_alleles/code.parallel/bindToEnv.R")

## Start up a parallel cluster
parallelCluster <- parallel::makeCluster(4) # parallelCluster <- parallel::makeCluster(parallel::detectCores())
print(parallelCluster)

print("<allele_counting>")

## savePropSchedule dictates that the results (in a list) will be saved every X*100% of the total number of sites
## i.e. 40000 with savePropSchedule = 0.05 will save results into an Rdata file every 2000 sites
## This is necessary for sequence data with large number of sites such as EBV, because the iterative process of addressing all possible pairs fails due to memory allocation, also becoming slower as more sites are added to the results

## Cut the data using an entropy minimum threshold
entropyNorm <- function(x){
    counts <- table(x)
    probs <- data.frame(AAs = names(counts), 
                        prob = as.numeric(counts/sum(counts)), 
                        stringsAsFactors=FALSE)
    en<--sum( probs$prob * log(probs$prob))
    return(en)
    }

entropies <- apply(data, 2, entropyNorm) # Apply entropyNorm function to dataset
removeIn <- which(entropies < entropyTrheshold) # Apply entropy cutoff 

if(length(removeIn)>0){ # Remove sites
    data= data[,-removeIn]
    dataNaming<- originaldataNaming[-removeIn]
    }else{
        data = data
        dataNaming <- originaldataNaming
    }

print(paste("<removed>",length(removeIn),"sites with entropy <=",entropyTrheshold)) 

dimension <- dim(data)
nCols <- dimension[2]
nRows <- dimension[1]
print(paste("<using data dimensions>","rows",nRows,"cols",nCols))


## Functions that will eventually be running in parallel
allelesCount <- function(x){
    alleles<- apply(cbind(data[,x[1]], data[,x[2]]), 1, paste, collapse = lociSeparator)
    return(list(table(alleles)))
    }
allelesNaming<- function(x){
    allelesNames<- paste(dataNaming[x[1]], dataNaming[x[2]], sep='-')
    return(allelesNames)
    }

print("<now allele counting...>")

out_allelesList <- list()
aa <- c()
saveEvents <- 0


# For each column/site get all possible pairs (i,j) with i=that column and j all columns >i, then try to parcel these pairs into groups and run the code in parallel
if(!is.null(nCols)){
    if(nCols > 0){
        pb <- txtProgressBar(min = 0, max = ncol(data)-1, style = 3)
        for(cc in 1:(ncol(data)-1)){
            allPairsColsT <- rbind(rep(cc,length((cc+1):ncol(data))),(cc+1):ncol(data))
            if((ncol(data)-cc) >= 5){ # If enough to parallelize
                ### Test 
                ## parcels <- round(seq(1,11,length.out=5))
                ## ranges <- list(parcels[1]:parcels[2], (parcels[2]+1):parcels[3], (parcels[3]+1):parcels[4], (parcels[4]+1):parcels[5])
                parcels <- round(seq(1, ncol(allPairsColsT), length.out=5)) # Parcels
                allPairsCols1<- as.matrix(allPairsColsT[,parcels[1]:parcels[2]], nrow=2) 
                allPairsCols2<- as.matrix(allPairsColsT[,(parcels[2]+1):parcels[3]], nrow=2) 
                allPairsCols3<- as.matrix(allPairsColsT[,(parcels[3]+1):parcels[4]], nrow=2)
                allPairsCols4<- as.matrix(allPairsColsT[,(parcels[4]+1):parcels[5]], nrow=2)
                allPairs<- list(allPairsCols1,allPairsCols2,allPairsCols3,allPairsCols4)
                parl<-TRUE
                }else{
                    allPairs<- list(allPairsColsT) # Not enough, put all pairs together
                    parl<-FALSE
                    }
                # debugNames <- function(x){ (dataNaming[x]) }
                # aa = rbind(aa, t(apply(allPairsColsT, 2, debugNames)))
                # Here, make sure to name all variables and functions that will be used per 'thread'
            mkWorker <- function() {
                bindToEnv(objNames=c('allelesNaming','allelesCount','reference','data','dataNaming','lociSeparator'))
                function(pairsParcel) {
                    pairsParcel = as.matrix(pairsParcel)
                    allelesList <- apply(pairsParcel, MARG = 2, FUN = allelesCount)
                    allelesList <- do.call(list, unlist(allelesList, recursive=FALSE))
                    allelesNames <- apply(pairsParcel, 2, allelesNaming)
                    names(allelesList)<- allelesNames
                    return(allelesList)
                  }
                }
            
            ## Run the code in parallel, in which the different parcels of pairs will be looked into to get their allele combinations
            n_sol_allelesList <- parallel::parLapply(parallelCluster, X=allPairs, fun=mkWorker())
            
            ## Save results depending on how many they are, given a parallel or not parallel run
                if(parl){
                  out_allelesList <- c(out_allelesList, n_sol_allelesList[[1]], n_sol_allelesList[[2]], n_sol_allelesList[[3]], n_sol_allelesList[[4]])
                }else{
                  out_allelesList<- c(out_allelesList, n_sol_allelesList[[1]])
                }
                if(cc %% (1+round(savePropSchedule*(ncol(data))))==0) {
                  saveEvents<- saveEvents+1  # Save this part of the results now
                  filename<-paste0(allelesOutFolder,tag,".N",saveEvents,".Rdata") 
                  save(out_allelesList, file=filename)
                  out_allelesList<- list()
                }
                setTxtProgressBar(pb, cc)
          }
          close(pb)
  }}


if(length(out_allelesList)!=0){
    saveEvents<- saveEvents+1 # There are still results in the variable not saved, save them now
    filename<-paste0(allelesOutFolder,tag,".N",saveEvents,".Rdata")
    save(out_allelesList, file=filename)
  }

## Save some useful variables for later
filename<-paste0(allelesOutFolder,tag,".extraInfo.Rdata")
save(saveEvents, file=filename)

## Clean results
out_allelesList<- list()

## Shutdown cluster neatly
  if(!is.null(parallelCluster)) {
    parallel::stopCluster(parallelCluster)
    parallelCluster <- c()
  }
