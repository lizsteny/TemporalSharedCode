#### Libraries
require(Biostrings)
require(seqRFLP)

################## Fasta to dataframe
fastaSeq <- readAAStringSet("denv4/data/seqs/denv4_aa.fasta")
id = names(fastaSeq)
seq = paste(fastaSeq)
denvseq <- data.frame(id, seq)

### Read csv file
denvtable <- read.csv("denv4/data/tables/denv4_table.csv")

### Merge
denvtable <- merge(denvseq,denvtable, by="id", all=F)
denvtable<-denvtable[-4]
denvtable<-data.frame(denvtable$id,denvtable$seq,denvtable$collection_date,denvtable$country,denvtable$name)

### Change names
names(denvtable) <- c("id", "seq", "collection_date", "country", "description")

### Subset
genotypeI<-denvtable[denvtable$genotype %in% denvtable$genotype[denvtable$genotype=="I"],]
genotypeII<-denvtable[denvtable$genotype %in% denvtable$genotype[denvtable$genotype=="II"],]

################## From csv to fasta
### Genotype I
genotypeI<-genotypeI[1:2]
genotypeI.fas = dataframe2fas(genotypeI, file="/Users/lizvillabona/Documents/Projects/denv_project/denv4/curation_data/denv4_genotypeI.fasta")

### Genotype II
genotypeII<-genotypeII[1:2]
genotypeII.fas = dataframe2fas(genotypeII, file="/Users/lizvillabona/Documents/Projects/denv_project/denv4/curation_data/denv4_genotypeII.fasta")




### Libraries
require(dplyr)
require(ggplot2)

################## Epitope table with frequency 
### Read epitopes csv
epitopes <- read.csv2("denv4/curation_data/tables/epitope_raw.csv")

### We only need starting and ending position of each epitope
epitopes <- epitopes[5:6]

### Create an empty list of vec
list_of_vecs<-list()

### Loop to calculate total number of epitopes at each position
for (i in 1:nrow(epitopes)){
    y <- c(epitopes$Starting.Position[i]:epitopes$Ending.Position[i]) # If you have in the 1st row a position with "5:8" in <y> you will have a vector with "5,6,7,8" 
    vecname <- paste('vec_', i, sep = '')
    assign(vecname,y) # Assign vecname to each vec in "y" 
    list_of_vecs[[i]]<-get(vecname)
}

### Calculate epitopes freq
epitopes<-bind_rows(lapply(list_of_vecs, as.data.frame.list))
epitopes <- t(epitopes)
epitopes<-table(epitopes)
epitopes <- as.data.frame(epitopes)
names(epitopes) <- c("locus", "freq")

region <- rep(c("Capsid", "Membrane" , "Envelope" , "NS1", "NS2A", "NS2B", "NS3", "NS4A", "NS4B", "NS5"),times=c(205,74,495,352,218,130,618,150,245,902))
locus <- c(1:3389)
denv <- cbind(region,locus) ; rm(locus)
epitopes <- merge(epitopes, denv, by =locus, all = FALSE)

### Save it
#write.csv(epitopes,"denv4/data/tables/denv4_loci_freq_IEDB.csv")


### epitopes MI
epitopes<-read.csv2("denv4/data/tables/denv4_loci_freq_IEDB.csv", stringsAsFactors = FALSE)
epitopes$freq

### Plot
ggplot(epitopes, aes(y = freq,
                     x = locus, 
                     color = genome.region)) +
    geom_point() +
    labs(y = "Frequency", 
         x = "Genome positions") +
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

#Capsid<-c(1:205)
#Membrane<-c(206:279)
#Envelope<-c(280:774)
#NS1<-c(775:1126)
#NS2A<-c(1127:1344)
#NS2B<-c(1345:1474)
#NS3<-c(1475:2092)
#NS4A<-c(2093:2242)
#NS4B<-c(2243:2487)
#NS5<-c(2488:3389)

### Read counts per H(x)
entropy$pieloueven<-as.numeric(as.character(entropy$pieloueven))

entropy$color <- factor(entropy$color, levels = c("Invariant", "Low", "Medium","High"))

entropy <- entropy %>%
    count(pieloueven,gen,color)

ggplot(entropy, aes(x = gen, 
                    y = n, 
                    fill = color)) + 
    geom_bar(position=position_fill(), 
             stat="identity") +
    theme_ipsum() +
    xlab("") +
    ylab("")  + 
    coord_flip() + 
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(name = "H(x)", values=c("#bdc9e1",
                                              "#74a9cf",
                                              "#2b8cbe",
                                              "#045a8d")) 

### Read csv
assays <- read.csv2("denv4/data/tables/denv4_epitope_assays.csv")
names(assays)

Bcells<-assays[ which(assays$assay=="B Cells Assays"),]
Tcells<-assays[ which(assays$assay=="T Cells Assays"),]
MHC<-assays[ which(assays$assay=="MHC ligand Assays"),]

ggplot(assays, aes(x = locus, 
                   y = hits, 
                   color = genome.region)) +
    xlab("Genome positions") + 
    ylab("Number of repetitions") + 
    geom_point(color="#95CF96", 
               shape=18, 
               alpha=0.5)  +
    theme_ipsum() 


### Read csv
barrita <- read.csv2("denv4/data/tables/denv4_loci_freq.csv")

barrita$hits.IEDB<-as.numeric(barrita$hits.IEDB)

### cuts
ggplot(barrita, aes(x = locus, y =hits.dataset, color=hits.IEDB))+
    xlab("Genome positions") + ylab("Number of repetitions") + 
    geom_point(alpha=0.9,size= 5,aes(colour = cut(Rep.IEDB, c(-1,0,1,2,5,10,15,20,25))))  + 
    theme_ipsum() +
    scale_color_manual(name = "Rep.IEDB", values=c("grey",
                                                   "#ccebc5",
                                                   "#a8ddb5",
                                                   "#7bccc4",
                                                   "#4eb3d3",
                                                   "#2b8cbe",
                                                   "#08589e"))

ggplot(barrita, aes(x = hits.dataset, y = hits.IEDB, color = genome.region)) + 
    labs(x="Scaled mutual information", y = "Linkage desequilibrium") + 
    geom_point(size = 5, 
               alpha = 0.9)  +
    geom_abline(intercept = 0.0, colour = "gray")  +
    theme_ipsum() +
    scale_color_manual(name = "Genome region", 
                       values=c( "#AADE94",
                                 "#fee08b", 
                                 "#8dd3c7", 
                                 "#80b1d3",
                                 "#b3de69", 
                                 "#bebada", 
                                 "#66c2a5")) +
    ylim(0,25) +
    xlim(0,8)
