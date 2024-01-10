########################################## Load libraries  ########################################## 
library(tidyverse)
library(lme4)
library(lmerTest) #("lmerTest") to get better pval from ANOVA-like test for Random Effects
library(DESeq2)
library(yaml) ## to load the config.yaml file
#setwd("/mnt/BioAdHoc/Groups/KronenbergLab/vcastelan/innateness/mouse/atac-seq/counts")
config <- read_yaml("config.yaml") ## reads yaml as a list
if(version$major == 4 & version$minor > 0.3) {
  setwd(config$innate_dir$server)
  outdir <- config$innate_output$server
} else {
  setwd(config$innate_dir$local)
  outdir <- config$innate_output$local
}
cat("working directory is: ", getwd(), "\n")
#Read count table with all cells
all_cells <- read.csv("data/countmatrix/ImmGen_ATACseq_count_table.csv", header = T, check.names = F) ## modify?

homer_index <- all_cells[,c(2,3,4,1)]
homer_index$not_used <- rep("NA",dim(homer_index)[1])
homer_index$strand <- rep("+",dim(homer_index)[1])
write.table(homer_index, file = "homer_index.bed", sep = "\t")

head(all_cells)
#assign peakID as rownames
rownames(all_cells) <- all_cells$PeakID


#Read metadata of the cells we are interested in
mdata <- read.csv("data/countmatrix/ImmGen_ATACseq_metadata.csv", header = T, check.names = F)
mdata <- mdata[-c(17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38),]
rownames(mdata) <- mdata$ID
mdata$ID <- NULL
#Get from count table just the cells we want to use (the ones in metadata)
#all_cells[,mdata$ID %in% colnames(all_cells)]
counts <- all_cells[,colnames(all_cells) %in% rownames(mdata)]
#Reorder columns of counts table to be in the same order as mdata
colnames(counts) 
rownames(mdata)
counts <- counts[,c(rownames(mdata))]

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = mdata,
                              design = ~ cell_type)

dds <- DESeq(dds)
res <- results(dds)
rld <- rlog(dds) # to get the log transform of gene expression
vsd <- varianceStabilizingTransformation(dds)
ntd <- normTransform(dds) #log2 transformation log2(n + 1)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$cell_type
colnames(sampleDistMatrix) <- NULL
ntop <- 1000 #number of MVG to take
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vsd)[select, ])
pca<-prcomp(mat,scale= TRUE)
pc_coord <- as.data.frame(pca$x[,])
data_ct <- as.data.frame(pca$x[,])
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))

#plot PCA
p <- ggplot(data_ct, aes(-pca$x[,1], -pca$x[,2], color=as.factor(mdata$cell_type))) +
  geom_point(size=5) +
  ggtitle("PC1 vs PC2") +
  labs(color = "Cell type")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw() + theme(text = element_text(size = 20))
ggsave('output/figures/pca_plot_atac.pdf')
##########################################  LMMs ########################################## 
mdata <- mdata[rownames(mdata) %in% row.names(pc_coord),] ## filter mdata for samples in pc_coord
gene_exp <- counts

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

master_table <- mdata
master_table$PC <- pc_coord$PC1
master_table <- cbind(master_table, t_gene_exp)
backup_master_table <- master_table


#Delete the genes without possible linear mixed model calculation
delete_genes<-c()
#for loop has to start where genes start in master table
for (i in 4:length(master_table)){
  suma<-sum(master_table[i])
  print(suma)
  if(suma == 0) { #Quit all genes with low expression because it will not be possible to calculate a lmm
    delete_genes<-rbind(delete_genes,i)
  }
}
delete_genes<-as.vector(delete_genes)
print("Genes that are going to be deleted:")
delete_gene_names<-colnames(master_table[delete_genes])
master_table[delete_gene_names] <- NULL
dg <- c("gene110754"
,"gene133147"
,"gene197873"
,"gene217242"
,"gene254046"
,"gene294882"
,"gene322096"
,"gene351612"
,"gene447353"
,"gene454980"
,"gene543805"
,"gene553271"
,"gene557827"
,"gene575884"
,"gene605305"
,"gene187616"
,"gene190448"
,"gene382397"
,"gene482089"
,"gene551977" 
,"gene469989"
,"gene527009"
,"gene582606"
,"gene105397"
,"gene167438"
,"gene169507"
,"gene175432"
,"gene215340"
,"gene237729"
,"gene302038"
,"gene327137"
,"gene392269"
,"gene459570"
,"gene501324"
,"gene541514"
,"gene590992"
,"gene596192")
master_table[dg] <- NULL

## create a model for each possible gene
all_genes <- colnames(master_table)[341649:length(master_table)]
print("############################### genes:")
print(length(master_table))
print("###############################qwertycode")
############# CALCULATE LMM ############# 
invisible(beta_lmm <- lapply(all_genes, function(k) {
  print(k)
  f = formula(paste(k,"PC + (1|cell_type)", sep = "~"))
  m = lmer(formula = f, data = master_table)
  m
}))

names(beta_lmm) <- all_genes # rm(beta_lmm)
############# SAVE R beta_lmm OBJECT ################ 
saveRDS(beta_lmm, file="beta_lmms3.rds")

beta_table <- data.frame(
  N = 1:length(all_genes),
  gene = all_genes,
  beta = sapply(beta_lmm, function(x){
    tryCatch(summary(x)$coefficients[2],
             error=function(e){NA})}), ## beta
  pval = sapply(beta_lmm, function(x){
    tryCatch(ranova(x)[["Pr(>Chisq)"]][[2]],
             error=function(e){NA})}), ## pval beta,
  varexp = sapply(beta_lmm, function(x){
    tryCatch(as.data.frame(VarCorr(x))$vcov[1] / (as.data.frame(VarCorr(x))$vcov[1] + as.data.frame(VarCorr(x))$vcov[2]),
             error=function(e){NA})}) ## variance explained
)
#write.csv(beta_table, file = "beta_table_pre_names3.csv")
### Generate a table guide with gene indexes and gene names equivalences
gene_guide <- data.frame(gene_symbol = gene_info$pseudo_name, gene_name = original_gene_names)
print("gene_guide")
gene_guide
rownames(gene_guide) <- gene_guide$gene_symbol
t_gene_guide <- as.data.frame(t(gene_guide))
t_gene_guide[delete_gene_names] <- NULL
gene_guide <- as.data.frame(t(t_gene_guide))
print("No trouble")

length(beta_table$gene)
length(gene_guide$gene_name)
tail(beta_table$gene)
tail(gene_guide$gene_name)

beta_table$gene <- gene_guide$gene_name
write.table(beta_table, file = "output/b_levels/results_mouse_beta_atac.tsv", sep = "\t")
