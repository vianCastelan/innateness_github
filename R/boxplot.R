#!/usr/bin/env Rscript
library(tidyverse)
#==============================# get data #==============================#
count <- read.csv(file = "data/countmatrix/immgen_ULI_RNAseq.csv", header = TRUE, check.names = FALSE)
colnames(count)[1] <- "gene"
# get metadata
metadata <- read.csv(file = "data/metadata/immgen_ULI_RNAseq_metadata.csv", header = TRUE)
rownames(metadata) <- metadata$names
metadata$X <- NULL
# filter counts and col_data
col_data <- metadata[metadata$tissue == "Spleen" & metadata$activation == "Naive", ]
counts <- count[, c("gene", col_data$names)]
#ranks for celltype
col_data$rank <- c(rep(6,5), rep(5,3), rep(1,6),rep(2,4), rep(4,2), rep(3,2))
#### plot per gene
gene_to_plot = as.character(commandArgs(TRUE)[1])
## plot 
merge(col_data,
      counts %>%
        filter(gene == gene_to_plot) %>%
        gather(key = "names", value = gene_to_plot)) %>%
  ggplot() +
  aes(x = reorder(all_conditions, rank),
      y = as.numeric(gene_to_plot), fill = all_conditions) +
  geom_boxplot() +
  theme_bw() + theme(text = element_text(size = 8), legend.position = "none") +
  ylab(gene_to_plot) + xlab(" ") + ggtitle(paste(gene_to_plot))
ggsave(filename = paste0("output/figures/boxplot_",gene_to_plot,".pdf"), width = 3, height = 3)
cat("file saved in output/figures")
# GENE="Nfat1c"
# Rscript R/boxplot.R $GENE