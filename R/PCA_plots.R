########## PCA plots
## Code from Vianka Cedillo Castelan

#===========# RNA-seq #===========#


## libraries
library(tidyverse)
library(yaml)
library(ggrepel)
library(DESeq2)
 ## reads yaml as a list 
config <- read_yaml("config.yaml")
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
#Do PCA analysis for spleen cells
colnames(counts) <- NULL
## DESeq2
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
ntop <- 1000 #number of MVG to take
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(assay(vsd)[select, ])
pca <- prcomp(mat,scale= TRUE)
#PCA plot
data_ct<-as.data.frame(pca$x[,])
percentVar <- round(100 * pca$sdev^2/sum(pca$sdev^2))

ggplot(data_ct, aes(-pca$x[,1], pca$x[,2], color=as.factor(col_data$cell_type))) +
  geom_point(size=6) +
  ggtitle("PC1 vs PC2") +
  labs(color = "Cell type")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() + theme(text = element_text(size = 20))

ggsave(filename = "output/figures/pca_plot.pdf", height = 6, width = 7)
ggsave(filename = "output/figures/pca_plot.png", height = 6, width = 7, dpi = 300)

#PC1 scores in each cell type
ggplot(data_ct, aes(x = reorder(col_data$cell_type,-pca$x[,1]), y = -pca$x[,1], fill=as.factor(col_data$cell_type)))+
  #scale_y_reverse()+
  geom_boxplot() +
  ggtitle("PC1") +
  xlab(paste0("Cell type")) +
  ylab(paste0("PC1")) +
  #coord_flip()+
  theme_bw() +
  theme(legend.position = "none") + theme(text = element_text(size = 18))

ggsave(filename = "output/figures/pc1_plot.pdf", height = 3, width = 7)
ggsave(filename = "output/figures/pc1_plot.png", height = 3, width = 7, dpi = 300)

#===========# ATAC-seq #===========#

#setwd("/mnt/BioAdHoc/Groups/KronenbergLab/vcastelan/innateness/mouse/atac-seq/counts")
#Read count table with all cells
all_cells <- read.csv("data/countmatrix/ImmGen_ATACseq_count_table.csv", header = T, check.names = F)

homer_index <- all_cells[,c(2,3,4,1)]
homer_index$not_used <- rep("NA",dim(homer_index)[1])
homer_index$strand <- rep("+",dim(homer_index)[1])
write.table(homer_index, file = "output/data/homer_index.bed", sep = "\t")

head(all_cells)
#assign peakID as rownames
rownames(all_cells) <- all_cells$PeakID


#Read metadata of the cells we are interested in
mdata <- read.csv("data/metadata/ImmGen_ATACseq_metadata.csv", header = T, check.names = F)
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

dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)

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
ggplot(data_ct, aes(-pca$x[,1], -pca$x[,2], color=as.factor(mdata$cell_type))) +
  geom_point(size=5) +
  ggtitle("PC1 vs PC2") +
  labs(color = "Cell type")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw() + theme(text = element_text(size = 20))

ggsave("output/figures/pca_atac.pdf", height = 6, width = 7) ## run me plz!!
ggsave("output/figures/pca_atac.png", height = 6, width = 7, dpi = 300) ## run me plz!!

#PC1 scores in each cell type
ggplot(data_ct, aes(x = reorder(mdata$cell_type,-pca$x[,1]), y = -pca$x[,1], fill=as.factor(mdata$cell_type)))+
  #scale_y_reverse()+
  geom_boxplot() +
  ggtitle("PC1") +
  xlab(paste0("Cell type")) +
  ylab(paste0("PC1")) +
  #coord_flip()+
  theme_bw() +
  theme(legend.position = "none") + theme(text = element_text(size = 18))

ggsave(filename = "output/figures/pc1_atac.pdf", height = 3, width = 7)
ggsave(filename = "output/figures/pc1_atac.png", height = 3, width = 7, dpi = 300)
