#!/usr/bin/env Rscript

## libraries
# library(ggvenn)
# require(gridExtra)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
myCol <- brewer.pal(3, "Pastel2")
## shared datasets
human <- unname(unlist(read.table("data/human/human_innateness_genes.tsv", sep = "\t", header = TRUE)))
mouse <- read.table("output/b_levels/results_mouse_beta_table_filt.csv", sep = ",", header = TRUE)
#mouse <- mouse[mouse$pval < 0.05 &  abs(mouse$beta) > 0.2,]
mouse <- mouse$gene
gene_comp <- list(mouse = mouse, human = human)
## venn object
# my_venn <- ggvenn(gene_comp,
#   fill_color = c("#0073C2FF", "#EFC000FF"),
#   stroke_size = 0.5, set_name_size = 4)
# grid.arrange(my_venn, top="Genes between human and mouse innateness gradiants", bottom="1000 MVG")
# ggsave("output/figures/venn_human_mouse.pdf")

## venn
myCol <- myCol[-3]
venn.diagram(
  x = gene_comp,
  category.names = c("Mouse", "Human"),
  filename = NULL, #'output/figures/venn_human_mouse.svg',
  output = TRUE,
  imagetype="png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  #lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #rotation = 1
) -> vnobj
ggsave(vnobj, file = "output/figures/venn_mouse_human.pdf")


## shared datasets
rna <- read.table("output/b_levels/results_mouse_beta_table_filt.csv", sep = ",", header = TRUE)
#rna <- rna[rna$pval < 0.05 & abs(rna$beta) > 0.2,]
atac <- read.table("output/b_levels/results_mouse_beta_atac.tsv", sep = "\t", header = TRUE)
rna <- rna$gene
atac <- atac$gene
gene_comp <- list(transcripts = rna, open_chromatin = atac)
## venn object
# my_venn <- ggvenn(gene_comp,
#   fill_color = c("#0073C2FF", "#EFC000FF"),
#   stroke_size = 0.5, set_name_size = 4)
# grid.arrange(my_venn, top="Genes shared between ATAC-seq and RNA-seq generated innateness models", bottom="1000 MVG")
# ggsave("output/figures/venn_rna_atac_all.pdf")
myCol <- brewer.pal(3, "Pastel2")
myCol <- myCol[-1]
venn.diagram(
  x = gene_comp,
  category.names = c("RNA", "Open Chromatin"),
  filename = NULL,
  output = TRUE,
  imagetype="png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  #lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #rotation = 1
) -> vnobj
ggsave(vnobj, file = "output/figures/venn_rna_atac.pdf")
