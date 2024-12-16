### Libraries
require(ggtree)
require(ggplot2)
require(ape)
require(hrbrthemes)
require(phytools)
require(dplyr)
require(stringr)

### Read GeomPS output
denvGeomPS <- read.csv2("denv4/data/output_geomps/geomps_rejecting001Xs.csv", stringsAsFactors = F)
denvGeomPS$siteA <- gsub("p","", denvGeomPS$siteA)
denvGeomPS$siteB <- gsub("p","", denvGeomPS$siteB)

denvMInorm <- read.csv("denv4/data/output_MI/TESTrun_Xs0.4_Al0_En0.4MI.csv", stringsAsFactors = F)

denvtable <- merge(denvGeomPS, denvMInorm, by = c("siteA","siteB"),all=F) ; rm(denvMInorm, denvGeomPS)
denvtable$MInorm <- as.numeric(as.character(denvtable$MInorm))
denvtable$R2 <- as.numeric(as.character(denvtable$R2))
denvtable$PSGEOMAB <- as.numeric(as.character(denvtable$PSGEOMAB))
denvtable$Sscore <- c(denvtable$MInorm*denvtable$PSGEOMAB)

names(denvtable) <- c("siteA","siteB", "PolA", "PSA", "PolB", "PSB",   "GeoMeanPS", "MI", "normMI", "morphA",  "morphB", "wMI", "wMInorm","wmorphA" ,"wmorphB" , "Dprime",  "R2", "Sscore")

denvtable$morph[is.na(denvtable$R2)] <- "Polymorphic"
denvtable$morph[is.na(denvtable$morph)] <- "Dimorphic"

### Density plot of normMI	
ggplot(denvtable, aes(x = normMI)) + 	
    geom_density() + 	
    theme_minimal() + 	
    labs(x = "NormMI")	

### Counting plot of normMI	
ggplot(denvtable, aes(x = normMI)) +
    geom_histogram(color = "white",
                   size = 0.15, aes(fill = morph)) +
    theme_ipsum()  +
    labs(x = "normMI") +
    scale_fill_manual(name = "",
                      values=c( "#A6BDD9",
                                "#8A68B0")) 
### normMI vs R2
ggplot(denvtable, aes(x = normMI, 
                      y = R2)) +
    labs(x="normMI", 
         y = "R2") +
    geom_point(size = 4,
               alpha = 0.7,
               colour="#A6BDD9")  +
    geom_abline(intercept = 0.0, colour = "grey") +
    expand_limits( x= c(0,1), y = c(0,1)) +
    theme_ipsum()

### Cor
cor.test(denvtable$normMI,denvtable$R2, method = "pearson")
quantile(denvtable$normMI)

### normMI vs geomps
ggplot(denvtable, aes(x = normMI, 
                      y = GeoMeanPS, 
                      color = Sscore)) +
    labs(x = "normMI", 
         y = "GeoMeanPS") +
    geom_point(size = 3.5,
               alpha = 0.9,
               aes(colour = cut(Sscore, c(-1,2,4,6,8,10,12,14,16)))) +
    theme_ipsum()  +
    scale_color_manual(name = "Sscore",
                       values=c("#ece2f0",
                                "#d0d1e6",
                                "#a6bddb",
                                "#3690c0",
                                "#02818a",
                                "#016c59",
                                "#014636")) +
    ylim(0,16)

### Density plot of GeoMeanPS	
ggplot(denvtable, aes(x = GeoMeanPS)) + 	
    geom_density() + 	
    theme_minimal() + 	
    labs(x = "GeoMeanPS")	

### Counting plot of GeoMeanPS	
ggplot(denvtable, aes(x = GeoMeanPS)) +
    geom_histogram(color = "white",
                   size = 0.15, aes(fill = morph)) +
    theme_ipsum()  +
    labs(x = "GeoMeanPS") +
    scale_fill_manual(name = "",
                      values=c( "#A6BDD9",
                                "#8A68B0"))

### Look at quantiles
quantile(denvtable$GeoMeanPS)
summary(denvtable$GeoMeanPS)

### Density plot of Sscore	
ggplot(denvtable, aes(x = Sscore)) + 	
    geom_density() + 	
    theme_minimal() + 	
    labs(x = "Sscore")	

### Counting plot of Sscore	
ggplot(denvtable, aes(x = Sscore)) +
    geom_histogram(color = "white",
                   size = 0.15, aes(fill = morph)) +
    theme_ipsum()  +
    labs(x = "Sscore") +
    scale_fill_manual(name = "",
                      values=c( "#A6BDD9",
                                "#8A68B0"))

### Look at the quantiles
quantile(denvtable$Sscore)
summary(denvtable$Sscore)

### Filter for S score
denvtable<-denvtable[denvtable$Sscore %in% denvtable$Sscore[denvtable$Sscore > 8],]

### Save final positions
#write.csv2(denvtable, "denv4/data/tables/denv4_final_loci_Xs001.csv", row.names = F)

### Read tree and table with the information
my.tree <- read.tree("denv4/data/output_ML/denv4.tree")
my.tree <- midpoint.root(my.tree)
table <- read.csv2("denv4/data/tables/denv4_table.csv",stringsAsFactors = FALSE)

### Get ids from tips
tips.names <- my.tree$tip.label
## Order info table
orden <- data.frame(tips.names)
orden$orden <- c(1:length(tips.names))
ldata <- length(tips.names)
names(orden) <- c("id","orden")
table <- merge(orden,table,by="id")
table <- arrange(table, orden)
table <- table[-2]

### Get info from tips
location <- table$country
names(location) <- table$id
region <- table$region
names(region) <- table$id
dates <- table$collection_date
names(dates) <- table$id
genotype <- table$genotype
names(genotype) <- table$id
tips.info <- data.frame(tips.names,dates,location,genotype)
tips.info <- as_tibble(tips.info)

### Create the vectors for the genotypes 
## fix to auto
I <- table[ which(table$genotype=='I'),]
II <- table[ which(table$genotype=='II'),]

I <- I$id
II <- II$id

gen <- list(I,II)
names(gen)<-c("I","II")

### Create vectors for region
America <- table[ which(table$region=='America'),]
AsiaOceAfri <- table[ which(table$region=='Asia+Oceania+Africa'),]

America <- America$id
AsiaOceAfri <- AsiaOceAfri$id

region <- list(America,AsiaOceAfri)
names(region)<-c("America","AsiaOceAfri")

### Plotting tree
mytree <- groupOTU(my.tree, region)

ggtree(mytree, layout = "circular", aes(color = group)) %<+% tips.info +
    scale_colour_manual(name = "Location",
                        values = c("#8da0cb","#a6d854"))

mytree <- groupOTU(my.tree, gen)

ggtree(mytree, layout = "circular", aes(color = group)) %<+% tips.info +
    scale_colour_manual(name = "Genotype",
                        values = c("#8da0cb","#a6d854"))

################## Trees

### Plotting the combinations of pairs of aminoacid
denv<-table[1:2]
my.number.cols<- length(unlist(strsplit(denv$seq[1], "")))

### Create names
my.col.names <- paste("", 1:my.number.cols, sep="p") 
denv <- str_split_fixed(string=denv$seq, pattern="", n=my.number.cols)
denv <- data.frame(id=table$id, p=denv)
names(denv) <- c("id", my.col.names)

denvtable <- denvtable[order(denvtable$GeoMeanPS, denvtable$normMI),]

denvtable <- denvtable %>% mutate(GeoMeanPS = format(GeoMeanPS, digits=4))
denvtable <- denvtable %>% mutate(normMI = format(normMI, digits=3))

sitioA<-denvtable$siteA
sitioB<-denvtable$siteB
sitios<-data.frame(sitioA,sitioB)
sitios$sitioA<-paste0("p",sitios$sitioA)
sitios$sitioB<-paste0("p",sitios$sitioB)
GeoMeanPS<-denvtable$GeoMeanPS
GeoMeanPS<-paste0("GeoMeanPS: ",GeoMeanPS)
GeoMeanPS<-paste0(GeoMeanPS, " ")
normMI<-denvtable$normMI
normMI<-paste0("normMI: ",normMI)


### Replace amino acid combinations below a threshold with XX  
## For example, lets say that loci12 (locus1xlocus2) formed combinations AB=10, CB=15 and AD=3. After applying a cutoff =< 3, loci12 will have combinations AB=10, CB=15 and XX=3. 
thresholdX <- 0.1 # X!X cutoff

genotype.count<-list()
region.count <- list()

for (p in 1:nrow(sitios))  { # Row sites for the loop 
    pp <- sitios[p,] # n number of row site 
    
    for (p1 in as.character(pp[1])) {
        paa1 <- denv[p1] } # Position of aa 1 / locus 1
    for (p2 in as.character(pp[2])) {
        paa2 <- denv[p2] } # Position of aa 2 / locus 2
    
    paa <- data.frame(denv$id,paste0(paa1[,1], paa2[,1])) # Loci 12 (locus1xlocus2)
    names(paa)<-c("id", "p") # Assign names in loci comb 
    tab <- plyr::count(paa[,2]) # Count aa
    nom <- tab$x # Name of aa
    freq <- tab$freq # Freq of aa
    
    for (f in 1:length(nom)) { # f in length of total aa names
        if (freq[f] <= ldata*thresholdX) {   # Replace aa combinations below n to XX
            paa$p <- gsub(nom[f], "XX", paa$p) } 
    } 
    
    n <- levels(as.factor(paa$p))
    l <- length(n)
    output <- list()
    
    for (i in 1:l) {
        aa <- select(filter(paa, paa$p == n[i]),id)
        num <- table(paa$p)
        aa <- as.character(aa$id)
        output [[i]] <- aa } # plot output
    
    n=paste0(n," : ",num)
    names(output) <- c(paste0(n))
    
    nam <- c(pp$sitioA)
    nam <- paste(nam, pp$sitioB)
    nam <- paste0(nam, " ")
    
    
    denv_ta<-read.csv2("denv4/data/tables/denv4_table.csv", stringsAsFactors = F)
    denv_gen<-denv_ta[1:3]
    denv_gen<-denv_gen[-2]
    denv_gen<-merge(paa,denv_gen, by="id", all=F)
    
    table.countgen <- denv_gen %>%
        count(genotype,p) 
    
    table.countgen %>% mutate(prob = prop.table(n))
    
    table.countgen$tree <- nam  # Maybe you want to keep track of which iteration produced it?
    genotype.count[[p]] <- table.countgen # Add it to your list for genotype %
    
    denv_re<-denv_ta[1:5]
    denv_re<-denv_re[-2]
    denv_re<-merge(paa,denv_re, by="id", all=F)
    
    table.countre <- denv_re %>%
        count(region,p) 
    
    table.countre %>% mutate(prob = prop.table(n))
    
    table.countre$tree <- nam  # Maybe you want to keep track of which iteration produced it?
    region.count[[p]] <- table.countre
    
    
    pdf(paste("denv4/data/output_trees/tree_", nam, ".pdf", sep = ""),width = 10, height = 8)
    
    geomps<-GeoMeanPS[p]
    dat<-paste0(geomps,normMI[p])
    nam <- paste0(nam, dat)
    
    dendogram <- force.ultrametric(my.tree, method = "extend")
    dendogram <- as.hclust(dendogram)
    dendogram <- as.dendrogram(dendogram)
    
    clus <- as.vector(paa$p)
    names(clus) <- paa$id
    
    g <- split(names(clus), clus)
    
    pl <- ggtree(dendogram, 
                 linetype='solid', 
                 layout = "circular", 
                 size=0.2) +  
        scale_color_manual(values = c("#af8dc3","#b2df8a","#fc8d59"))
    
    #clades <- sapply(g, function(x) MRCA(pl, x)) they don't form clades
    
    pl <- groupOTU(pl, g, group_name = "AA") +
        aes(color=AA)
    
    rownames(table) <- table$id
    
    d <- data.frame(label = names(clus),
                    region = table[names(clus), "region",],
                    genotype = table[names(clus), "genotype",])
    row.names(d) <- d$label; d$label <- NULL
    
    tree <- gheatmap(pl, d,width=0.1,colnames=FALSE) + 
        scale_fill_manual(values=c("#e0f3f8", "#fee08b", "#d7191c", "#2c7bb6"), name="")+
        ylab(paste0(nam))
    
    print(tree)
    dev.off()
} 

genotype.count = do.call(rbind, genotype.count)	
genotype.count$prob <- genotype.count$n/ldata*100
region.count = do.call(rbind, region.count)	
region.count$prob <- region.count$n/ldata*100

#write.csv2(genotype.count, "denv4/data/output_trees/genotype_count.csv", row.names = F)
#write.csv2(region.count, "denv4/data/output_trees/region_count.csv", row.names = F)

## Clean GE
rm(list = ls(all.names = TRUE)) 
