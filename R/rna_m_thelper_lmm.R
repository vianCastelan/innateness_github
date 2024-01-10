#!/usr/bin/env Rscript

## Linear Mixed Models for Innateness Score of mouse Immune cells from ImmGen Datasets
## Gabriel Ascui & Vianka Cedillo Castellan
## 25JAN21  ///// MOD 14FEB21

########################################## Load libraries  ########################################## 
library(tidyverse)
library(lme4)
library(lmerTest) #### install.packages("lmerTest") to get better pval from ANOVA-like test for Random Effects
library(parallel)
##########################################  load data     ########################################## 
# col_data <- read.csv(file = "../../master_tpm_metadata_2020-12-01_justnaive_th1_th17ranking.csv",header=TRUE) ## still has activated T cells
# counts <- read.csv(file = "Th1_count_table.csv",header = TRUE, check.names = FALSE) ## DESeq2 by Vianka

count <- read.csv(file = "data/countmatrix/immgen_ULI_RNAseq.csv", header = TRUE, check.names = FALSE)
head(count)
rownames(count) <- count[, 1]
# get metadata
metadata <- read.csv(file = "data/metadata/immgen_ULI_RNAseq_metadata.csv", header = TRUE)
rownames(metadata) <- metadata$names
metadata$X <- NULL
metadata$names <- NULL
##
print("Data has been read")

#==============================# T helper and ILCs to create new LMM #==============================#
# filter counts and col_data
col_data <- na.omit(metadata)
counts <- count[, rownames(col_data)]
dim(counts)

########################################################################################################
#########################################   LMM pre processing   #######################################
########################################################################################################
# pseudo names to avoid problem with RIKEN genes on columns of the master table 
gene_info <- data.frame(gene_symbol = row.names(counts), pseudo_name = NA)
gene_info$pseudo_name <- sapply(100001:sum(100000,nrow(counts)), function(x) paste0("gene",x))
######## master table with all genes 
## gene names
temporal_gene_names <- gene_info$gene_symbol
row.names(counts) <- make.names(gene_info$gene_symbol, unique=TRUE)
#counts <- counts[,-1]
## replace gene names by pseudo names 
row.names(counts) <- gene_info$pseudo_name
t_counts <- t(counts) 
master_table <- col_data
master_table <- cbind(master_table, t_counts)
# gene = colnames(master_table)[7] ## Genes start at position 7
#Quit all genes with gen exp=1
delete_genes <- c()
for (i in 9:length(master_table)){
  suma<-sum(master_table[i])
  #print(suma)
  if(suma==0) {
    delete_genes <- rbind(delete_genes,i)
  }
}
delete_genes <- as.vector(delete_genes)
cat("Genes that are going to be deleted:","\n")
delete_gene_names <- colnames(master_table[delete_genes])
cat(delete_gene_names, "\n")
master_table[delete_gene_names] <- NULL
## put model in a list for iterations
all_genes <- colnames(master_table)[8:length(master_table)]


### make linear mixed models function
lmm_model <- function(g = all_genes, tab, xvar = "Th1") {
    require(lme4)
    require(lmerTest)
    require(parallel)
    print("initiating LMM calculations")
    ############# THIS STEP TAKES AROUND 40 MIN IN MY PERSONAL COMPUTER #############
    invisible(model_list <- mclapply(g, function(k) {
        #print(k) #
        f <- formula(paste(k, paste0(xvar, " + (1|cell_type)"), sep = "~"))
        m <- lmer(formula = f, data = tab)
    }, mc.cores = 16))
    names(model_list) <- g # rm(model_list)
    ############# SAVE R MODEL_LIST OBJECT ################
    print("Models are done!")
    saveRDS(model_list, file = paste0("output/data/", xvar, "_lmm.rds"))
    print("Models are saved!")
    beta_table <- data.frame(
        N = seq_along(all_genes),
        gene = all_genes,
        beta = sapply(model_list, function(x) {
            tryCatch(summary(x)$coefficients[2],
                error = function(e){NA})}), ## beta
        pval = sapply(model_list, function(x) {
            tryCatch(ranova(x)[["Pr(>Chisq)"]][[2]],
                error = function(e) {NA})}), ## pval beta,
        varexp = sapply(model_list, function(x) {
            tryCatch(as.data.frame(VarCorr(x))$vcov[1] / (as.data.frame(VarCorr(x))$vcov[1] + as.data.frame(VarCorr(x))$vcov[2]),
                error = function(e) {NA})}) ## variance explained
    )
    ###### WRITE CSV WITH BETA TABLE ####
    ### Generate a table guide with gene indexes and gene names equivalences
    gene_guide <- data.frame(gene_symbol = gene_info$pseudo_name, gene_name=temporal_gene_names)
    rownames(gene_guide) <- gene_guide$gene_symbol
    t_gene_guide <- as.data.frame(t(gene_guide))
    t_gene_guide[delete_gene_names] <- NULL
    gene_guide <- as.data.frame(t(t_gene_guide))
    beta_table$gene <- gene_guide$gene_name
    write.table(beta_table, file = paste0("output/b_levels/", xvar,"_beta_table.tsv"), sep = "\t")
}

#### Run function for each pathway:
lmm_model(g = all_genes, tab = master_table, xvar = "Th1")
lmm_model(g = all_genes, tab = master_table, xvar = "Th2")
lmm_model(g = all_genes, tab = master_table, xvar = "Th17")