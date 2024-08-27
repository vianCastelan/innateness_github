#!/usr/bin/env Rscript

## libraries
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
## source
source("R/IPA_functions.R")
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

##=================================================== Run for loop
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

## Annotations for pheatmap
rld <- rld[ ,-c(1:9)] ## remove ILCs 
col_data <- col_data[-c(1:9), ] ## remove ILCs 

# heatmap function
zscore_heatmap <- function(pathway, path_name, filt_max = NA, s = FALSE) {
  require(pheatmap)
  require(RColorBrewer)
  ## Access rld and col_data from the global environment
  rld <- get("rld", envir = .GlobalEnv)
  col_data <- get("col_data", envir = .GlobalEnv)
  ## get config info
  if (any(grepl("rld",ls())) & any(grepl("col_data",ls()))) {
    ##----- annotations and colors
    ann <- data.frame(celltype = factor(c(rep(6,5),rep(4,3),rep(1,6),rep(2,4),rep(5,2),rep(3,2)), labels = c("CD4","CD8","MAIT","NKT","gdT","NK")))
    rownames(ann) <- colnames(rld)
    Groups <- read.csv("colors.csv")$color
    names(Groups) <- as.vector(read.csv("colors.csv")$category)
    ann_colors <- list(celltype = Groups[1:6])
    ##----- process mat 
    mat <- rld[rownames(rld) %in% pathway, ] - rowMeans(rld[rownames(rld) %in% pathway, ]) # rld filter and get z-score
    if (!is.na(filt_max)) {
      mat <- mat[apply(mat, 1, function(x) max(x)) > filt_max, ]
    }
    if (nrow(mat) == 0) {
      print("Matrix is empty after filtering...")
      return(NULL)
    }
    mat <- mat[ , order(ann$celltype)]
    ##----- row labels 
    row_labs <- lapply(rownames(mat), function(x) bquote(italic(.(x))))
    ##----- pheatmap function
    pheatmap(
      mat,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
      clustering_callback = callback,
      cluster_cols = FALSE,
      annotation_col = ann,
      annotation_colors = ann_colors,
      cellheight = 10, cellwidth = 10,
      filename = if (s == TRUE) {
        paste0("output/figures/", path_name, "_heatmap.pdf")
      } else { 
        NA
      },
      width = 6,
      height = 6, 
      labels_col = NA,
      labels_row = as.expression(row_labs)
    ) -> p
    return(p)
  } else {
    print("rld or col_data not found in environment")
    return(NULL)
  }
}

## apply function
## pathways 
paths <- list(
  NK_signaling = c(
    "Fcer1g",
    "Ncr1",
    "Nck1",
    "Klra7",
    "Klrc1",
    "Klrb1",
    "Klrd1",
    "Syk"
  ),
  Integrins = c(
    "Itgal",
    "Itgb1",
    "Itgb2",
    "Itgax",
    "Itga4",
    "Itgam",
    "Itgae"
  ),
  Cytokine_receptors = c(
    "Il2rb",
    "Il12rb2",
    "Il18rap",
    "Il18r1",
    "Stat4",
    "Stat3",
    "Il10ra",
    "Il10rb"
  ),
  PI3K = c(
    "Pik3ca",
    "Pik3cb",
    "Pik3cd",
    "Pik3cg",
    "Pik3r1",
    "Pik3r5",
    "Akt1",
    "Akt2",
    "AKt3",
    "Foxo1",
    "Foxo3",
    "Plcg1",
    "Plcg2"
  ),
  ERK = c(
    "Mapk3",
    "Mapk1",
    "Map2k1",
    "Map2k2",
    "Map3k1",
    "Ets1",
    "Elk1",
    "Jun",
    "Junb",
    "Jund",
    "Fos",
    "Fosb",
    "Fosl1",
    "Fosl2"
  ),
  TCR = c(
    "Icam1",
    "Lfa1",
    "Cd28",
    "Lat",
    "Lck",
    "Itk",
    "Lcp2",
    "Fyn",
    "Nfatc1",
    "Ppp3cb",
    "Ppp3ca",
    "Zap70",
    "Camk4",
    "Camk2g",
    "Calm1"
  ),
  EIF2 = c(
    "Rps2",
    "Rps20",
    "Eif5",
    "Eif4a2",
    "Eif3k",
    "Eif3e",
    "Rpl4",
    "Rpl8"
  )
)

#
zscore_heatmap(pathway = paths$NK_signaling, path_name = "NK_signaling", s = TRUE)
zscore_heatmap(pathway = paths$Integrins, path_name = "Integrins_signaling", s = TRUE)
zscore_heatmap(pathway = paths$Cytokine_receptors, path_name = "Cytokine_receptors", s = TRUE)
zscore_heatmap(pathway = paths$PI3K, path_name = "PI3K_signaling", s = TRUE)
zscore_heatmap(pathway = paths$ERK, path_name = "ERK_signaling", s = TRUE)
zscore_heatmap(pathway = paths$TCR, path_name = "TCR_signaling", s = TRUE)
zscore_heatmap(pathway = paths$EIF2, path_name = "EIF2_signaling", s = TRUE)

## pathways
can <- read.table("output/ipa/canonical_pathways_all.txt", sep = "\t", skip = 2, header = TRUE)
##--- get first pathway (function from IPA_functions.R)
# Natural Killer Cell signaling
pathway <- process_string(can[2,5]) # 2
# TCR signaling pathway
pathway <- process_string(can[3,5]) # 3
# Motility
pathway <- process_string(can[20,5]) # 20
# ERK signaling
pathway <- process_string(can[41,5]) # 41
# PI3K pathway
pathway <- process_string(can[23,5]) # 23
# PKCtheta signaling pathway
pathway <- process_string(can[22,5]) # 22
# EIF2 signaling pathway
pathway <- process_string(can[1,5]) # 1


