#!/usr/bin/env Rscript

## libraries
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
#==============================# get data #==============================#
count <- read.csv(file = "data/countmatrix/immgen_ULI_RNAseq.csv", header = TRUE, check.names = FALSE)
head(count)
rownames(count) <- count[, 1]
count <-  count[, -1]
# get metadata
metadata <- read.csv(file = "data/metadata/immgen_ULI_RNAseq_metadata.csv", header = TRUE)
rownames(metadata) <- metadata$names
metadata$X <- NULL
metadata$names <- NULL
# filter counts and col_data
col_data <- metadata[metadata$tissue %in% c("Spleen","Small_Intestine") & metadata$activation == "Naive", ]
counts <- count[, rownames(col_data)]

## regular log transform data log10(x + 1 )
rld <- DESeq2::rlog(as.matrix(counts)) ## this takes a while
rownames(rld) <- rownames(counts)
############# pheatmaps function
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

### Reference paper
load("output/data/human_Thelper.RData")

## Run for loop
genes <- list(th1 = th1_genes,th2 = th2_genes, th17 = th17_genes)
for (i in seq_along(genes)){
    mat <- rld[genes[[i]],] - rowMeans(rld[genes[[i]],]) ### this converts it into fold_change, how much does it move from the mean of the row.
    pheatmap(mat,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
    clustering_callback = callback,
    cutree_rows = 5,
    filename = paste0("output/figures/", names(genes[i]), "_heatmap.pdf"),
    width = 7,
    height = 18)
}
