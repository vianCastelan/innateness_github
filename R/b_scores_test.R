#!/usr/bin/env Rscript

## load libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
## load data
b_rna <- read.table("output/b_scores/immgen_scores.tsv", sep = "\t", header = TRUE) ## no "#" characters
b_atac <- read.table("output/b_scores/innateness_scores_atac2.tsv", sep = "\t", header = TRUE)
b_human <- read.table("output/b_scores/innateness_scores_human.tsv", sep = "\t", header = TRUE)
b_mait <- read.table("output/b_scores/ssMAITcells_sc_scores.tsv", sep="\t", header = TRUE)
b_salmonella <- read.table("output/b_scores/salmonella_sc_scores.tsv", sep="\t", header = TRUE)



## load metadata
meta <- read.csv("data/metadata/immgen_ULI_RNAseq_metadata.csv") #
meta_salmonella <- read.csv("data/metadata/metadata_salmonella.csv")
meta_mait <- read.csv("data/metadata/steady_state_metadata.csv")

## remove additional activated CD4
meta <- meta[-c(56,57),]

## combine data
comb <- data.frame(
    names = meta$names,
    cell_type = meta$cell_type,
    mouse_b = b_rna$innateness_score[-c(53,54,55)],
    human_b = b_human$innateness_score[-c(53,54,55)],
    atac_b = b_atac$innateness_score[-c(53,54,55)],
    activation = meta$activation,
    tissue = meta$tissue
)

## plot: Mouse vs. Human models
comb %>%
    filter(tissue == "Spleen", activation == "Naive") %>%
    ggplot() +
    aes(x = mouse_b, y = human_b) +
    geom_smooth(method = "lm",
                se = TRUE,
                col = "grey",
                fill = "light grey",
                alpha = .5) +
    geom_point(aes(col = cell_type, size = 3)) +
    #geom_text(aes(label = names)) +
    theme_bw() #+ scale_x_log10() + xlim(-100000000, 100000000)
ggsave("output/figures/comp_h_m.pdf", height = 5, width = 6)

## plot: RNA vs. ATAC models
comb %>%
    filter(tissue == "Spleen", activation == "Naive") %>%
    ggplot() +
    aes(x = mouse_b, y = -atac_b) +
    geom_smooth(method = "lm",
                se = TRUE,
                col = "grey",
                fill = "light grey",
                alpha = .5) +
    geom_point(aes(col = cell_type, size = 3)) +
    #geom_text(aes(label = names)) +
    theme_bw() #+ 
    #ylim(0,5000)
    #scale_x_log10() + xlim(-100000000, 100000000) +
    #scale_y_log10() + ylim(0, 3000)
ggsave("output/figures/comp_rna_atac.pdf", height = 5, width = 6)


### many tests: make a long graph

## load innateness score data data
isaac <- read.table("output/b_scores/GSE157225_scores.tsv", sep = "\t", header = TRUE)
flavell <- read.table("output/b_scores/flavell01_scores.tsv", sep = "\t", header = TRUE)
goldrath <- read.table("output/b_scores/goldrath01_scores.tsv", sep = "\t", header = TRUE)
## plot NKT
options(scipen=10000)
isaac %>%
    separate(cell_type, into = c("celltype","tissue"), sep = ":") %>%
    mutate(tissue = str_remove(tissue, "\\..*")) %>%
    filter(tissue == "spleen", celltype != "CD4_T", celltype != "NKT") %>% 
    ggplot() +
    aes(x = innateness_score, 
        y = reorder(celltype, innateness_score), 
        fill = reorder(celltype, innateness_score)) + 
    geom_boxplot() + theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) -> NKT

## flavell
flavell %>% 
    ggplot() + 
    aes(x = innateness_score,
        y = reorder(cell_type, innateness_score),
        fill = reorder(cell_type, innateness_score)) +
    geom_point(size=5) + theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) -> Th17

## goldrath
meta <- read.csv("data/metadata/goldrath_metadata.csv")
goldrath %>%
    merge(.,meta, by.x="cell_type", by.y = "cell") %>% 
    filter(!cell_type %in% c("GSM235554.CEL","GSM235555.CEL")) %>% ## weird points 
    ggplot() +
    aes(x = innateness_score,
        y = reorder(cell_type.y, innateness_score),
        fill = reorder(cell_type.y, innateness_score)) +
    geom_boxplot() + theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) + scale_x_log10() + xlim(0,300)-> CD8

## ssMAIT cells
b_mait$cell_type <- gsub("\\.","-",b_mait$cell_type)
mait_t <- merge(b_mait, meta_mait, by.x = "cell_type", by.y = "X")
mait_t %>% filter(seurat_clusters %in% c(0,1,3,5,6,10)) %>%
    ggplot() +
    aes(x = reorder(seurat_clusters, innateness_score), 
        y = innateness_score, fill = factor(seurat_clusters)) +
    geom_violin() + 
    geom_boxplot(width = 0.1, fill = "white") +
    theme_bw() + coord_flip() + 
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14)) -> ss_MAIT



## Salmonella MAIT cells
library(Seurat)
load("../mait_libs/salmonella.RData")
DimPlot(salmonella, reduction = "tsne", group.by = "orig.time_point") 
DimPlot(salmonella, reduction = "tsne", label = TRUE)
## cluster 1 Klrg1+ 
## cluster 3,4,9 Il7r+


b_salmonella$cell_type <- gsub("\\.","-",b_salmonella$cell_type)
mait_s <- merge(b_salmonella, meta_salmonella, by.x = "cell_type", by.y = "X")
mait_s %>% filter(orig.time_point == "DAY40") %>% 
    filter(seurat_clusters %in% c(1,3,4,9)) %>%
    mutate(day_40_cat = ifelse(seurat_clusters %in% c(3,4,9),"Il7r+","Klrg1+")) %>%
    ggplot() + 
    aes(x = reorder(day_40_cat, -innateness_score), 
        y = innateness_score) +
    geom_violin() + 
    geom_boxplot(width = 0.05, fill = "white") + 
    theme_bw() +
    theme(axis.text = element_text(size = 14)) -> innate_salmonella


## including salmonella day0 and day6
mait_s %>% filter(orig.time_point %in% c("DAY40","DAY6","Uninfected")) %>%
    filter((seurat_clusters == 13 & orig.time_point == "Uninfected") | (seurat_clusters %in% c(10,7,8) & orig.time_point == "DAY6") | (seurat_clusters %in% c(1,3,4,9) & orig.time_point == "DAY40")) %>%
    ggplot() +
    aes(x = reorder(orig.time_point, innateness_score),
        y = innateness_score) +
    geom_violin() + 
    geom_boxplot(width = 0.05, fill = "white") + 
    theme_bw() +
    theme(axis.text = element_text(size = 14)) -> innate_salmonella2

## plot grid 
plot_grid(Th17,CD8,NKT,ss_MAIT,innate_salmonella,innate_salmonella2, ncol = 3)
ggsave("output/figures/b_scores_applied.pdf", width = 18, height = 8)

