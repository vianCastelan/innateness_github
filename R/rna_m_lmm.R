########################################## Load libraries  ########################################## 
cat("loading R libraries","\n")
library(tidyverse)
library(lme4)
library(lmerTest) #("lmerTest") to get better pval from ANOVA-like test for Random Effects
library(DESeq2)
library(yaml) ## to load the config.yaml file
########################################## read configuration  ########################################## 
config <- read_yaml("config.yaml") ## reads yaml as a list
if(grepl("/mnt/",getwd())) {
  setwd(config$innate_dir$server)
  outdir <- config$innate_output$server
} else {
  setwd(config$innate_dir$local)
  outdir <- config$innate_output$local
}
cat("working directory is: ", getwd(), "\n")
##################################### read arguments from bash  ####################################### 
## Rscript R/rna_h_lmm.R $pval $beta
filt_pval = as.numeric(commandArgs(TRUE)[1])
filt_beta = as.numeric(commandArgs(TRUE)[2])
########################################## load data & metadata ########################################## 
## Read data from "data" directory
count <- read.csv(file = "data/countmatrix/immgen_ULI_RNAseq.csv",header = TRUE, check.names = FALSE)
rownames(count) <- count[,1]
count$gene_symbol <- count[1]
cat("dimension of count matrix are: ",dim(count), "\n")
## read metadata
metadata <- read.csv(file = "data/metadata/immgen_ULI_RNAseq_metadata.csv", header=TRUE)
rownames(metadata) <- metadata$names
metadata$X <- NULL
metadata$names <- NULL
#Create the variables that are going to be used in the models
col_data <- metadata[metadata$tissue == "Spleen" & metadata$activation == "Naive",]
gene_exp <- count[c("gene_symbol", rownames(col_data))]
counts <- gene_exp
counts$gene_symbol <- NULL
colnames(counts) <- NULL
cat("filtered dimensions for non-activated splenocytes are: ", dim(counts))
cat("Data has been read and prepared for downstream analysis. Initiating DESeq2...","\n")
########################################## DESeq2 ########################################## 
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = col_data,
                              design = ~ cell_type)
dds <- DESeq(dds)
res <- results(dds)
rld <- rlog(dds) # to get the log transform of gene expression
vsd <- varianceStabilizingTransformation(dds)
ntd <- normTransform(dds) #log2 transformation log2(n + 1)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$all_conditions
colnames(sampleDistMatrix) <- NULL
ntop <- config$n_mvg #number of MVG to take (1000 as default)
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vsd)[select, ])
pca<-prcomp(mat,scale= TRUE)
pc_coord <- as.data.frame(pca$x[,])
cat("PC coordinates have been generated...","\n")
##########################################  filter samples ########################################## 
col_data <- col_data[rownames(col_data) %in% row.names(pc_coord),] ## filter col_data for samples in pc_coord
gene_exp %>% select(gene_symbol, rownames(col_data)) -> gene_exp ## filter gen_exp for samples in col_data
gene_exp[1] <- NULL
########################################################################################################
########################################  every gene LMM loop ##########################################
######################################################################################################## 
cat("Initiating Linear Mixed Modeling for each gene.", "\n")
# pseudo names to avoid problem with RIKEN genes on columns of the master table 
gene_info <- data.frame(gene_symbol = row.names(gene_exp), pseudo_name = NA)
gene_info$pseudo_name <- sapply(100001:sum(100000,nrow(gene_exp)), function(x) paste0("gene",x))
######## master table with all genes 
## gene names
original_gene_names <- rownames(gene_exp)
## replace gene names by pseudo names 
row.names(gene_exp) <- gene_info$pseudo_name
#create a master table with all the info (tissue, cell type, activation, and gene expression), cell names in rows and genes in columns
t_gene_exp <- t(gene_exp)
master_table <- col_data
master_table$PC <- pc_coord$PC1
master_table <- cbind(master_table, t_gene_exp)
#Delete the genes without possible linear mixed model calculation
delete_genes<-c()
#for loop has to start where genes start in master table
for (i in 9:length(master_table)){
  suma<-sum(master_table[i])
  #print(suma)
  if(suma < dim(master_table)[1]) { #Quit all genes with low expression because it will not be possible to calculate a lmm
    # This consideres that a gene is expressed if the sum is higher than the number of samples, 
    # as if the expression is 1 or less in all samples.
    delete_genes<-rbind(delete_genes,i)
  }
}
delete_genes<-as.vector(delete_genes)
cat("Genes that are going to be deleted:",length(delete_genes),"\n")
delete_gene_names<-colnames(master_table[delete_genes])
master_table[delete_gene_names]<-NULL
## create a model for each possible gene
all_genes <- colnames(master_table)[9:length(master_table)]
############# CALCULATE LMM ############# 
cat("running linear models...","\n")
invisible(beta_lmm <- lapply(all_genes, function(k) {
  f = formula(paste(k,"PC + (1|cell_type)", sep = "~"))
  m = lmer(formula = f, data = master_table)
  m
}))
names(beta_lmm) <- all_genes # rm(beta_lmm)
############# SAVE R beta_lmm OBJECT ################ 
saveRDS(beta_lmm, file=paste0(outdir,"/b_levels/","beta_lmm.rds"))
##### Create data frame with all the derived data from beta_lmm
beta_table <- data.frame(
  gene = names(beta_lmm),
  beta = sapply(beta_lmm, function(x){
    tryCatch(summary(x)$coefficients[2]*-1,
             error=function(e){NA})}), ## beta (slope)
  sd = sapply(beta_lmm, function(x){
    tryCatch(sqrt(diag(vcov(x)))[2],
             error=function(e){NA})}), ## standard deviation of beta
  pval = sapply(beta_lmm, function (x){
    tryCatch(summary(x)$coefficients[10],
             error=function(e){NA})}), ## pval for beta from model
  pval_anova = sapply(beta_lmm, function(x){
    tryCatch(anova(x)[["Pr(>F)"]],
             error=function(e){NA})}), ## pval beta from anova test
  cov = sapply(beta_lmm, function(x){
    tryCatch(vcov(x)[2,2],
             error=function(e){NA})}),
  varexp = sapply(beta_lmm, function(x){
    tryCatch(as.data.frame(VarCorr(x))$vcov[1] / (as.data.frame(VarCorr(x))$vcov[1] + as.data.frame(VarCorr(x))$vcov[2]),
             error=function(e){NA})}), ## variance explained
  singular = sapply(beta_lmm, function(x){
    tryCatch(isSingular(x),
             error=function(e){NA})}) ## possible type I error
)
## normalized beta values : 
library(scales)
beta_table$beta_norm <- rescale(beta_table$beta, to=c(-2,2))
### WRITE CSV WITH BETA TABLE ####)
### Generate a table guide with gene indexes and gene names equivalences
gene_guide <- data.frame(gene_symbol = gene_info$pseudo_name, gene_name = original_gene_names)
rownames(gene_guide) <- gene_guide$gene_symbol
t_gene_guide <- as.data.frame(t(gene_guide))
t_gene_guide[delete_gene_names] <- NULL
gene_guide <- as.data.frame(t(t_gene_guide))
if (length(beta_table$gene) == length(gene_guide$gene_name)) print("No trouble") else print("ups... there is a problem...")
tail(beta_table$gene)
tail(gene_guide$gene_name)
## write beta table 
beta_table$gene <- gene_guide$gene_name
write.table(beta_table, file = paste0(outdir,"/b_levels/results_mouse_beta_table.tsv"), sep = "\t")
write.csv(beta_table, file = paste0(outdir,"/b_levels/results_mouse_beta_table.csv"))
## filtered table 
write.csv(beta_table %>% filter(pval < 0.05), file = paste0(outdir, "/b_levels/results_mouse_beta_table_filt.csv"))
cat("RNA linear mixed models and beta table are all done!","\n")
