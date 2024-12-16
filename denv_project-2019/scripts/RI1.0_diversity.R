######### Diversity index plots

### Libraries
library(stringr)
require(vegan)
require(ggplot2)
require(hrbrthemes)
require(Cairo)

### Read csv file
# for denv1 change to "denv1/"
denvtable <- read.csv2("denv4/data/tables/denv4_table.csv", stringsAsFactors = FALSE)

### Only ID and genomic sequence
denv <- denvtable[1:2]

### Names cols
my.number.cols<- length(unlist(strsplit(denv$seq[1], "")))

### Create names
my.col.names <- paste("", 1:my.number.cols, sep="") 

### Separate positions
denv <- str_split_fixed(string=denv$seq, pattern="", n=my.number.cols)
denv <- data.frame(id=denvtable$id, p=denv)
names(denv) <- c("id", my.col.names)
rm(denvtable)

### Filtered seqs
alleles<-denv[-1]

### Normalized shannon index
pielouEvennessCalc<- function(X){
  sam<- alleles[,which(colnames(alleles)==X)]
  gaps<- which(!(sam %in% c('-','X')))
  sam<- sam[gaps]
  Y<- t(as.matrix(table(sam)))
  n<- length(unique(sam))
  d<- diversity(Y, index = "shannon")/log(n)
  return(d)
}

### Shannon index
shannonCalc<- function(X){
  Y<- t(as.matrix(table(alleles[,which(colnames(alleles)==X)])))
  d<- diversity(Y, index = "shannon")
  return(d)
}

### Number of polymorphism
polysCalc<- function(X){
  sam<- alleles[,which(colnames(alleles)==X)]
  n<- length(unique(sam))
  return(n)
}

shan<- unlist(lapply(colnames(alleles), FUN=shannonCalc))
poly<- unlist(lapply(colnames(alleles), FUN=polysCalc))
peve<- unlist(lapply(colnames(alleles), FUN=pielouEvennessCalc))

### Create data frame with the diversity index
denv<- data.frame(ID=colnames(alleles), shannon=shan, npoly=poly, pieloueven=peve, stringsAsFactors=FALSE)

### Change names
names(denv) <- c("id", "shan", "poly","pieloueven")

### Nan = 0
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
denv[is.nan(denv)] <- 0

denv$id <- as.numeric(denv$id)
denv$shan <- as.numeric(denv$shan)
denv$pieloueven <- as.numeric(denv$pieloueven)
denv$poly <- as.numeric(denv$poly)

denv$genome.region <- rep(c("Capsid", "Membrane" , "Envelope" , "NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5"),times=c(205,74,495,352,218,130,618,150,245,902))

### Open pdf to save graphics
#CairoPDF("denv4/graphics/denv1_diversity.pdf", width = 12, height = 7)
# for denv1 change to "denv1/"
# for only one genotype change the .pdf

### Plot number of polymorphisms per position
ggplot(denv, aes(x = id,
                 y = poly, 
                 color = genome.region)) +
  labs(y = "Number of polymorphisms", 
       x = "Genome positions",
       title = "Number of polymorphisms per position for DENV4") +
  geom_point(alpha = 0.8) +
  theme_ipsum() +
    scale_color_manual(name = "Genome region", 
                       values = c("#AADE94",
                                "#BB92BE", 
                                "#fee08b", 
                                "#8dd3c7", 
                                "#80b1d3", 
                                "#fdb462", 
                                "#b3de69", 
                                "#bebada", 
                                "#3288bd", 
                                "#66c2a5")) 

### Plot shannon index
ggplot(denv, aes(x = id,
                 y = shan, 
                 color = genome.region)) +
  labs(y = "Values", 
       x="Genome positions", 
       title = "Shannon index for DENV4") +
  geom_point() +
  theme_ipsum() +
  scale_color_manual(name = "Genome region", 
                     values = c("#AADE94",
                                "#BB92BE", 
                                "#fee08b", 
                                "#8dd3c7", 
                                "#80b1d3", 
                                "#fdb462", 
                                "#b3de69", 
                                "#bebada", 
                                "#3288bd", 
                                "#66c2a5")) 

### Plot normalized shannon index
ggplot(denv, aes(x=id,
                 y=pieloueven, 
                 color = genome.region)) +
  labs(y = "Values", x="Genome positions", 
       title = "Normalized shannon index for DENV4") +
  geom_point(alpha = 0.8) +
  theme_ipsum() +
    scale_color_manual(name = "Genome region", 
                       values=c("#AADE94",
    "#BB92BE", 
    "#fee08b", 
    "#8dd3c7", 
    "#80b1d3", 
    "#fdb462", 
    "#b3de69", 
    "#bebada", 
    "#3288bd", 
    "#66c2a5")) 

#dev.off()

### Write csv with diversity information
# for denv1 change to "denv1/"
write.csv(denv,"denv4/data/output_diversity/denv4_diversity.csv", row.names = F)

## Clean GE
rm(list = ls(all.names = TRUE)) 
