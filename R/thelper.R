#!/usr/bin/env Rscript
library(biomaRt)
library(stringr)
#==============================#  read table_file for human T helper signatures #==============================# 
table_genes <- read.table("data/human/th_table_human.csv", header = T, check.names = F, sep = ",", skip = 1) ### find file
table(table_genes$CellType)
th1_ens <- table_genes$Gene[table_genes$CellType == "CD4_TH1"]
th2_ens <- table_genes$Gene[table_genes$CellType == "CD4_TH2"]
th17_ens <- table_genes$Gene[table_genes$CellType == "CD4_TH17"]
count <- read.csv(file = "data/countmatrix/immgen_ULI_RNAseq.csv", header = TRUE, check.names = FALSE)
rownames(count) <- count[, 1]
#==============================# BioMart #==============================#
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

IDs_th1 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
             filters = "ensembl_gene_id", values = th1_ens,
             mart = mart)
th1_genes <- str_to_title(unique(IDs_th1$hgnc_symbol))

IDs_th2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", values = th2_ens,
                 mart = mart)
th2_genes <- str_to_title(unique(IDs_th2$hgnc_symbol))

IDs_th17 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", values = th17_ens,
                 mart = mart)
th17_genes <- str_to_title(unique(IDs_th17$hgnc_symbol))


th1_genes <- th1_genes[th1_genes %in% rownames(count)]
th2_genes <- th2_genes[th2_genes %in% rownames(count)]
th17_genes <- th17_genes[th17_genes %in% rownames(count)]

save(th1_genes, th2_genes, th17_genes, file = "output/data/human_Thelper.RData")
