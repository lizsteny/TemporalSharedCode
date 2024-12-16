print("<allele_counts_to_tables>")

filename<-paste0(allelesOutFolder,tag,".extraInfo.Rdata")
load(filename) # Load saveEvents

if(saveEvents > 0){
      for(Nsave in 1:saveEvents){

              print(paste("<converting list to matrix-like> N",Nsave,"of",saveEvents))

              filename <- paste0(allelesOutFolder, tag, ".N", Nsave, ".Rdata")
              load(filename) # Load out_allelesList
              alleleNames <- names(out_allelesList)

{
listout <- list()
listout1 <- list()
{
    for (l in 1:length(out_allelesList)) {
        out_df<-as.data.frame(out_allelesList[l])
        comb<-as.character(out_df[,1]) # Look possible combinations
        all_count<-as.numeric(out_df[,2]) # Count of that possible combinations
        {
            for (i in 1:nrow(out_df))  { 
                if (all_count[i] <= threholdX) { # Apply threholdX cutoff
                    comb[i] <- "X!X" }  } # Replace in comb with X!X (change if other loci separator is used)
            out_df[1] <- comb # Replace with changes
            if (all_count[i] <= threholdX) {
                XX <- out_df[out_df[,1] == "X!X",]
                XXrow <- c("X!X", sum(as.numeric(XX[,2])))
                out_df <- out_df[out_df[,1]!= "X!X",]
                out_df <- rbind(out_df, XXrow) }
                out_df[,1] <- gsub("-", "X", out_df[,1]) # Gaps are assumed as any state X
                }
        listout<-append(listout,list(out_df))
        listout1<-append(listout1,list(as.vector(t(out_df))))
        }
    }
names(listout)<-alleleNames
names(listout1)<-alleleNames
}
              
out_allelesList <- listout
final_allelesTable <- listout1
final_alleleNames<- names(out_allelesList)

## Now some allele combinations may be XX because they were XX, or --, X-, -X
## Same can happen with others, naturally happening XP and -P will become both XP
## These possible repetitions of XX and XX or XP and XP need to me collapsed
print("<collapsing repetead combination due to - -> X>")
scanRep2Collapse <- function(X){
    allelesChar <- X[seq(1, length(X), 2)]
    allelesCharTab <- table(allelesChar)
    allelesRep <- names(allelesCharTab)[which(allelesCharTab>1)]
    for(arep in allelesRep){
        inds <- which(X == arep)
        if(length(inds > 1)){
        countsXX<- X[inds + 1]
        sumCountsXX<- sum(as.numeric(countsXX))
        X<- X[-c(inds, inds + 1)]
        X<- c(X, arep, sumCountsXX)
        }
    }
    return(X)
}

# Add progress bar and use function that sums X!X comb in one
coll_final_allelesTable<- pblapply(final_allelesTable, FUN=scanRep2Collapse) 

## Save the corrected table
final_allelesTable <- coll_final_allelesTable
print("<formating entries>")

## Build final table entry              
allLengths <- sapply(final_allelesTable, length)

## Compute maximum length of allele combinations found for each pair
maxLength <- max(allLengths)
              
## Print("<fill in empty spaces>") and add empty values values to list elements to have complete size matrix
allelesMatrix <- lapply(final_allelesTable, function(v) { c(v, rep('', maxLength-length(v)))})

## Print("<transform to matrix>") and force matrix on it
allelesMatrix <- do.call(rbind, allelesMatrix)

## Add the site pairs
alleleMatrixNames <- matrix(unlist(strsplit(final_alleleNames,'-')), ncol = 2, byrow = TRUE)
allelesMatrix<- cbind(alleleMatrixNames[,1], alleleMatrixNames[,2], allelesMatrix)

filename<- paste0(allelesOutFolder, tag, "_allele_table.N", Nsave, ".csv")
print(paste("<write to table> ", filename))
write.table(allelesMatrix, filename, col.names=FALSE, row.names=FALSE, sep=',')

## Clean some memory
allelesTable <- NULL
out_allelesList <- NULL
final_allelesTable <- NULL
final_alleleNames <- NULL
alleleMatrixNames <- NULL
allelesMatrix <- NULL
    }
}