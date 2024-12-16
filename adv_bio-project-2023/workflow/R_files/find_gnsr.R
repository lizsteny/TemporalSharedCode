#!/usr/bin/env Rscript
# Script: find_gnsr.R
# call: find_gnsr.R [path_tx2] [metadata] [logfc2] [FDR_thresh] [perc] \
#                   [filtered] [shared] [tableFC] [abundance]
# path_tx2: str, path to tx2 file for genome
# metadata: str, path to author metadata.txt
# logfc2: int, absolute log2fc cutoff
# FDR_thresh: float, FDR threshold
# filtered: str, path to output for filtered.tsv
# shared: str, path to output for shared.tsv
# tableFC: str, path to output for tableFC.tsv
# abundance: multiple strings, paths to all abundance.h5 files (separate by spaces)
# Authors: Liz Villabona-Arenas, Mathijs Balk, Rik Aikes 
# Purpose: Conducts differential expression analysis on gene expression data and filters results.
# Date: December, 2023

# make sure shared R libraries can be found
Path_shared_dir <- "/prog/BIF30806/project/groen_team/R_shared_libraries"
.libPaths(c(Path_shared_dir, .libPaths()))

# retrieve command line arguments when executing R script
args <- commandArgs(trailingOnly = TRUE)

# load libraries
library(tximportData)
library(tximport)
library(edgeR)
library(statmod)
library(limma)
library(rhdf5)

# create paths to abundance files for quantification results 
#h5.samples <- args[9,-1]
h5.samples <- args[-(1:8)]
# extract sample names from directory names
sample_names <- sub(".*/kallisto_A-(.*)/abundance.h5", "\\1", h5.samples)
# match with the experiment number
names(h5.samples) <- sample_names

# read from genome specific tx2gene
tx2gene <- read.table(args[1], header=TRUE, sep = '\t')

# for import transcripts-to-gene mapping and count data using tximport
txi.kallisto.tsv <- tximport(h5.samples, type = 'kallisto', tx2gene = tx2gene, ignoreAfterBar = TRUE)
countTable <- round(txi.kallisto.tsv$counts)

# import metadata from "metadata.txt" as in the paper
meta <- read.table(args[2],header=T)
meta <- meta[rep(1:nrow(meta),each=5),]

# make repeat column
meta$Repeat <- 1:5

# add R - repeats next to experiment number
meta[,1] <- paste(meta[,1],"R",1:5,sep="")

# make rownames based on col1 so that matches experiment names
rownames(meta) <- paste(meta[,1],sep="") 

# column names should match the samples in your count data
meta <- meta[colnames(countTable),]

# limma analysis
meta$Condition <- as.factor(meta$Condition) 
condition <- relevel(meta$Condition,ref="Axenic")
design <- model.matrix(~ condition)
colnames(design) <- sub("condition","",colnames(design))

# create a dgelist
dge <- DGEList(counts= countTable) 

# filtering low counts
good <- filterByExpr(dge, design)

# number of genes retained
dge <- dge[good,,keep.lib.sizes=FALSE]

# effective library sizes
# normalize counts using trimmed mean of M values (TMM) (adjust for differences in library sizes between samples)
dge <- calcNormFactors(dge)

# perform voom transformation
v <- voom(dge,design)

# fit linear model and perform eBayes moderation
fit <- lmFit(v, design)
fitEx <- eBayes(fit)

# extract information from fitEx with a focus on the contrast (coef)
# create an empty list to store filtered results
filtered_genes_list <- list()
logFC_threshold <- as.numeric(args[3])
FDR_threshold <- as.numeric(args[4])
perc_of_lists <- as.numeric(args[5])
# iterate over coefficients
for (i in 2:40) {
  # get the results for the current coefficient
  current_results <- topTable(fitEx, coef = i, number = Inf, sort.by = "P", adjust.method = "BH")
  
  # filter genes based on log-fold change (1) and FDR threshold (0.01)
  filtered_results <- current_results[abs(current_results$logFC) > logFC_threshold & current_results$adj.P.Val < FDR_threshold, ]
  
  # store the filtered results in the list
  filtered_genes_list[[i]] <- filtered_results
  
  # print or save logFC and FC for each coefficient (just to look)
  # cat("Coefficient", i, "\n")
  # cat("Number of genes:", nrow(filtered_results), "\n")
}

# extract row names from each dataframe in filtered_genes_list
row_names_list <- lapply(filtered_genes_list, rownames)

# get gene occurences across all dataframes
DEGs_occurrences <- table(unlist(row_names_list))
DEGs_occurrences <- as.data.frame(DEGs_occurrences)
names(DEGs_occurrences) <- c("DEGs","occurrences")

# filter genes that appear in at least 50% of the lists
filter_threshold <- 39 * perc_of_lists
selected_DEGs <- DEGs_occurrences$occurrence >= filter_threshold

# genes that appear more than the threshold
filtered_DEGs_occurrences <- DEGs_occurrences[selected_DEGs, ]
cat(nrow(filtered_DEGs_occurrences), "DEGs consistently show upregulation of at least",filter_threshold, "times across all treatments")
write.table(filtered_DEGs_occurrences, file = args[6], sep = "\t", quote = FALSE, row.names = FALSE)

# comparing to the paper
DEGs_paper <- c(
  "AT1G02920", "AT1G02930", "AT1G21110", "AT1G21120", "AT1G26380",
  "AT1G26410", "AT1G26420", "AT1G35230", "AT1G64170", "AT1G65500",
  "AT1G76930", "AT2G19190", "AT2G25470", "AT2G30750", "AT2G39200",
  "AT2G43620", "AT3G46280", "AT4G12490", "AT4G12500", "AT4G23140",
  "AT4G23220", "AT4G28420", "AT5G24110", "AT1G26390"
)

# check which genes from genes_of_paper are present in the row names of none_coef
shared_DEGs <- DEGs_paper[DEGs_paper %in% filtered_DEGs_occurrences[,1]]
cat(length(shared_DEGs), "DEGs have been identified as upregulated in both the mentioned paper and our analysis")
write.table(shared_DEGs, file = args[7], sep = "\t", quote = FALSE, row.names = FALSE)

# make vector with GNSR gene names (also convert factors back to strings)
GNSR_genes <- as.character(filtered_DEGs_occurrences[,1])

# table with log|FC|
tableFC <- topTable(fitEx, number = Inf, sort.by = "none", adjust.method = "BH")

# make selection for only GNSR genes
tableFC.GNSR <- tableFC[GNSR_genes,]
write.table(tableFC.GNSR, file = args[8], sep = "\t", quote = FALSE, row.names = TRUE)
