#!/usr/bin/env Rscript
# Script: graphs.R
# Authors: Geethanjali Srinivasan, Liz Villabona-Arenas, Mathijs Balk 
# Purpose: Generate heat maps and bar plots
# Date: December, 2023

##### heat maps

# load libraries
library(pheatmap)

# read GNSR file 
GNSR_Araport11 <- read.table("/prog/BIF30806/project/groen_team/results_deg/tableFC_GNSR_Araport11.tsv", header = T)
GNSR_TAIR10 <- read.table("/prog/BIF30806/project/groen_team/results_deg/tableFC_GNSR_TAIR10.tsv", header = T)

# subset 
GNSR_Araport11 <-GNSR_Araport11[,1:39]
GNSR_TAIR10 <- GNSR_TAIR10[,1:39]

# rownames association
group_info<-read.csv("group_info.csv", header = T)
rownames(group_info) <- group_info$X
group_info1 <- group_info[2]

# create a heatmap for Araport11
pheatmap(
  GNSR_Araport11,
  annotation_col = group_info1,
  cluster_rows = TRUE,
  scale = "row", 
  main = "Cluster analysis of DEGs - Araport11",
  color = colorRampPalette(c("#FFFFFF", "#C1DD87", "#7BC86D", "#39651F"))(100),
  angle_col = 90,  # rotate x-axis labels
  border_color = NA,
)

# create a heatmap for TAIR10
pheatmap(
  GNSR_TAIR10,
  annotation_col = group_info1,
  cluster_rows = TRUE,
  scale = "row", 
  main = "Cluster analysis of DEGs - Araport11",
  color = colorRampPalette(c("#FFFFFF", "#C1DD87", "#7BC86D", "#39651F"))(100),
  angle_col = 90,  # rotate x-axis labels
  border_color = NA,
)

##### bar plots
Counts_Araport11 <- read.table("/prog/BIF30806/project/groen_team/results_deg/Araport11_total_DEGs.tsv", header = T)
Counts_TAIR10 <- read.table("/prog/BIF30806/project/groen_team/results_deg/TAIR10_total_DEGs.tsv", header = T)

# load libraries
library(ggplot2)

# Combining the datasets
combined_data <- rbind(transform(Counts_Araport11, dataset = "Araport11"),
                       transform(Counts_TAIR10, dataset = "TAIR10"))

ggplot(combined_data, aes(x = reorder(condition, -DEGs), y = DEGs, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(
    x = "Treatments",
    y = "Number of DEGs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # rotate x-axis labels
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank(),   # remove minor grid lines
    axis.line = element_line(color = "black", size = 0.2),  # set axis line color to black and size to 0.2
    axis.ticks = element_line(color = "black", size = 0.2)  # set axis tick color to black and size to 0.2
  ) +
  guides(fill = guide_legend(title = NULL)) +  # remove legend title
  scale_fill_manual(values = c("Araport11" = "#66c2a5", "TAIR10" = "#AEDB78"))
