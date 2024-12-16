### Libraries
require(ggplot2)
require(plyr)
require(hrbrthemes)
require(Cairo)

## Read csv file
denvtable <- read.csv2("denv4/data/tables/denv4_table.csv", stringsAsFactors = FALSE)
# for denv1 change to "denv1/"

## Open pdf to save graphics
#CairoPDF("denv4/graphics/denv_general.pdf", width = 12, height = 8)

## Plot collection date of the genomic sequences per country by genotype
ggplot(denvtable, 
       aes(y=country,x=collection_date,colour=genotype)) +
  labs(x = "Collection year", y="Contries (n=23)", title="Collection date of the genomic sequences per country by genotype") +
  geom_point(alpha=0.6,size = 2.3) + 
  theme_ipsum() +
  scale_colour_manual(values = c("#8da0cb","#e78ac3","#ffd92f","#a6d854"))

## Create a dataframe only with genotype and collection_date
years<-denvtable[3:4]

# Count
years<-count(years)

# Plot number of genomic sequences per year
ggplot(years, 
       aes(collection_date,freq,color=genotype)) +
  labs(y = "Number of sequences", x="Collection date",title="Number of genomic sequences per year") + geom_point(alpha=0.6,size = 2.5) + 
  theme_ipsum() +
  scale_colour_manual(values = c("#8da0cb","#e78ac3","#ffd92f","#a6d854"))

#dev.off()


## Clean GE
#rm(list = ls(all.names = TRUE)) 
