### Libraries
require(dplyr)
require(phangorn)
require(ape)
require(phytools)
require(stringr)

### Read MI 
denvMInorm <- read.csv("denv4/data/output_MI/TESTrun_Xs0.4_Al0_En0.4MI.csv",stringsAsFactors = FALSE)
denvMInorm$MInorm <- as.numeric(as.character(denvMInorm$MInorm))

sitioA <- denvMInorm$siteA
sitioB <- denvMInorm$siteB
sitios <- data.frame(sitioA, sitioB)
sitios$sitioA <- paste0("p", sitios$sitioA)
sitios$sitioB <- paste0("p", sitios$sitioB)

### Read tree and table with the information
mstr.tree <- read.tree("denv4/data/output_ML/denv4.tree")

### Root tree from midpoint
mstr.tree <- midpoint.root(mstr.tree)
mstr.taxnames <- mstr.tree$tip.label

### Read denv table with seqs
table <- read.csv2("denv4/data/tables/denv4_table.csv",stringsAsFactors = FALSE)

### Plotting the combinations of pairs of aminoacids
denv <- table[1:2]
my.number.cols <- length(unlist(strsplit(denv$seq[1], "")))

### Create names
my.col.names <- paste("p", 1:my.number.cols, sep = "") 
denv <- str_split_fixed(string = denv$seq, pattern = "", n = my.number.cols)
denv <- data.frame(id = table$id, p = denv)
names(denv) <- c("id", my.col.names)
denv1 <- denv
rownames(denv1) <- denv1$id
denv1 <- denv1[-1]

AAmx <- as.matrix(denv1)
nAAs <- dim(sitios)[1]
totalSeq <- dim(AAmx)[1]
namemx <- rownames(denv1)

### Work out name mapping between sequences and tree
name.map <- vector("numeric",length=length(namemx))

for(i in 1:totalSeq){
    name.map[i] = which(namemx == mstr.taxnames[i])
}
rm(i)

### Reorder
reordered.AAmx <- denv1[name.map,]
reordered.names <- namemx[name.map]
rm(AAmx,namemx)

### Polymorphism determining functions
site.polymorphism <- function(aVector,threshold){
    temp <- table(aVector)
    return(length(which(temp >= threshold)))  
}

### For invariant sites
is.site.invariant <- function(aVector){
    temp <- table(aVector)
    nImpAAs = length(temp)
    if(nImpAAs == 1){
        return (1)
    } else{
        return (0)
    }
}

### Calculate empirical parsimony score per site
empirical.parsimony1 <- array(0,c(nAAs,3)) ; empirical.parsimony2 <- array(0,c(nAAs,3)) #changed to three

### Replace amino acid combinations below a threshold with XX  
## For example, lets say that loci12 (locus1xlocus2) formed combinations AB=10, CB=15 and AD=3. After applying a cutoff =< 3, loci12 will have combinations AB=10, CB=15 and XX=3.  
thresholdX <- 0.01 # X! cutoff
thresholdAlleles <- 0 # Reject alleles below this number per locus

{
    for (p in 1:nrow(sitios))  { # Row sites for the loop 
        pp <- sitios[p,] # n number of row site 
        {
            
        for (p1 in as.character(pp[1])) {
            paa1 <- reordered.AAmx[p1] } # Position of alleles in locus 1
            for (p2 in as.character(pp[2])) {
            paa2 <- reordered.AAmx[p2] } # Position of alleles in locus 2
            
            tab <- plyr::count(paa1[,1]) # Count alleles in locus 1
            nom <- tab$x # Name of alleles in locus 1
            freq <- tab$freq # Freq of alleles in locus 1
            
            tab2 <- plyr::count(paa2[,1]) # Count alleles in locus 2
            nom2 <- tab2$x # Name of alleles in locus 2
            freq2 <- tab2$freq # Freq of alleles in locus 2
            
            
            for (f in 1:length(nom)) { # f in length of total allele names
                if (freq[f] <= thresholdX*totalSeq) {  
                paa1[,1] <- gsub(nom[f], "X", paa1[,1]) } # Replace alleles below n to X
            }
            
            for (f2 in 1:length(nom2)) { # f in length of total aa names
                if (freq2[f2] <= thresholdX*totalSeq) {  
                    paa2[,1] <- gsub(nom2[f2], "X", paa2[,1]) } # Replace aa combinations below n to X
            }

        st1 <- as.vector(paa1[,1]) # Locus 1 with X alleles that formed XX comb
        names(st1) <- reordered.names
        st2 <- as.vector(paa2[,1]) # Locus 2 with X alleles that formed XX comb
        names(st2) <- reordered.names
        }
        
## Parsimony per locus!
        {
            if(is.site.invariant(st1)){ # If invariant site = 0
                P<-0 }
            else {
                stphydat <- phyDat(st1, type = "USER", levels = c("A","R","N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X"))
                names(stphydat) = reordered.names
                # Calculate parsimony steps for locus 1!
                P <- parsimony(mstr.tree,stphydat) }

            threshold <- thresholdAlleles*totalSeq # Reject rare alleles
            polymorphism <- site.polymorphism(st1,threshold)
            
            empirical.parsimony1[p,]<- c(p1,polymorphism,P) }
        {
            ##  Calculate parsimony steps for locus 2
            if(is.site.invariant(st2)){
                P<-0 }
            else {
                stphydat <- phyDat(st2,type="USER", levels = c("A","R","N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X"))
                names(stphydat) = reordered.names
                
                # Calculate the parsimony !
                P <- parsimony(mstr.tree,stphydat) }
            
            threshold <- thresholdAlleles # Reject rare alleles
            polymorphism <- site.polymorphism(st2,threshold)
            empirical.parsimony2[p,]<- c(p2,polymorphism,P) }
    } 
}
    
empirical.parsimony1 <- as.data.frame(empirical.parsimony1)
empirical.parsimony2 <- as.data.frame(empirical.parsimony2)
empirical.parsimony<-cbind(empirical.parsimony1,empirical.parsimony2)
colnames(empirical.parsimony) <- c('siteA','PolA','PSA','siteB','PolB','PSB')


empirical.parsimony$PSA <- as.numeric(as.character(empirical.parsimony$PSA))
empirical.parsimony$PSB <- as.numeric(as.character(empirical.parsimony$PSB))


### Geom mean cannot be calculated with 0 
gm_mean = function(x, na.rm = TRUE){
    exp(sum(log(x+1), na.rm = na.rm) / length(x))-1
}


### Calculate geometric mean of the parsimony score per loci pair
empirical.parsimony$PSGEOMAB <- apply(cbind(empirical.parsimony$PSA, 
                                            empirical.parsimony$PSB), 
                                      MARGIN = 1, 
                                      FUN = gm_mean)

### Write output
write.csv2(empirical.parsimony,"denv4/data/output_geomps/geomps_rejecting001Xs.csv",row.names = FALSE)

### Clean GE
rm(list = ls(all.names = TRUE)) 
