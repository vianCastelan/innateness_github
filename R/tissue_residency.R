#! usr/bin/env Rscript
#
# load libraries
library(Seurat) ## v5
library(ggplot2)
library(tidyverse)

if (!file.exists("data/GSE131847/exhaustion_join.RData")) {
  # load data 
  files <- list.files(path = "data/GSE131847/countmatrix", pattern = "*.tsv.gz", full.names = TRUE)
  sobj_list <- list()
  for (f in files) {
    fn <- gsub(".*/(.*?)\\.tsv\\.gz", "\\1", f)
    dtable <- read.table(f, header = TRUE, row.names = 1, check.names = FALSE)
    mat <- t(as.matrix(dtable))
    sobj_list[[fn]] <- CreateSeuratObject(counts = mat, assay = "RNA", project = fn)
    sobj_list[[fn]] <- NormalizeData(sobj_list[[fn]], verbose = FALSE)
    sobj_list[[fn]] <- FindVariableFeatures(sobj_list[[fn]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    
  }
  ## merge objects without integration
  #cells <- merge(x = sobj_list[[1]], y = sobj_list[2:length(sobj_list)])
  #cells
  # find anchors and integrate data
  anchors <- FindIntegrationAnchors(object.list = sobj_list)
  cells <- IntegrateData(anchorset = anchors)
  cells
  ## integration Seurat v5 
  # cells <- IntegrateLayers(cells, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")
  cells[["RNA"]] <- JoinLayers(cells[["RNA"]])
  ## add metadata from names
  cells@meta.data$orig.cells <- gsub(".*_([^_]+_[^_]+_cell).*", "\\1", cells@meta.data$orig.ident)
  cells@meta.data$orig.cells <- ifelse(cells@meta.data$orig.cells == "TRM_2_cell","D32_TRM_2_cell", cells@meta.data$orig.cells)
  cells@meta.data$orig.cells <- gsub(".*_([D]\\d+_cell)-.*", "\\1",cells@meta.data$orig.cells)
  cells@meta.data$orig.cells <- gsub(".*_(N\\w+_cell)-.*", "\\1",cells@meta.data$orig.cells)
  cells@meta.data$orig.tissue <- ifelse(grepl("TRM_cell|TRM_2_cell", cells@meta.data$orig.cells), "Gut","Spleen")
  cells@meta.data$orig.timepoint <- ifelse(grepl("Naive_cell", cells@meta.data$orig.cells), 
                                           0, 
                                           sub(".*D(\\d+)_.*", "\\1", cells@meta.data$orig.cells))
  sort(table(cells@meta.data$orig.timepoint))
  ## rm temp objects
  # rm(dtable,mat,sobj_list)
  ## post-processing
  # cells <- NormalizeData(cells)
  # cells <- FindVariableFeatures(cells)
  cells <- ScaleData(cells) ### this will scale integrated assay 
  cells@active.assay <- "RNA"
  ## save big tsv file and RMD file
  write.csv(x = as.matrix(GetAssayData(object = cells, layer = "counts", assay = "RNA")),
              file = "data/GSE131847/exhaustion.csv")
  save(cells, file = "data/GSE131847/exhaustion_int.RData")
  ## Check Variance Dimensionality
  cells <- RunPCA(cells)
  DimPlot(cells, group.by = "orig.tissue")
  # ElbowPlot(cells)
  # ## Clustering
  cells <- FindNeighbors(cells, dims = 1:10)
  cells <- FindClusters(cells, resolution = 0.5)
  # ## Dimensional Reduction
  cells <- RunUMAP(cells, dims = 1:10)
  DimPlot(cells, split.by = "orig.timepoint", label = TRUE)
  # ## apply innateness score
  b_exhausted <- read.table("output/b_scores/GSE141847_scores.tsv", sep="\t", header = TRUE)
  row.names(b_exhausted) <- b_exhausted$cell_type
  exhausted_s <- merge(b_exhausted, cells@meta.data, by = "row.names")
  write.csv(exhausted_s, "data/metadata/exhausted_metadata.csv")
  save(cells, file = "data/GSE131847/exhaustion_join.RData")
}

## start over
load("data/GSE131847/exhaustion_join.RData")
# meta <- read.csv("data/metadata/exhausted_metadata.csv", row.names = 2)
meta <- exhausted_s[-1,]
row.names(meta) <- meta$cell_type
cells@meta.data <- meta
# ## plot innateness
FeaturePlot(cells, features = "Klrg1")
FeaturePlot(cells, features = "innateness_score", split.by = "orig.timepoint")
VlnPlot(cells, features = "innateness_score", pt.size = 0, group.by = "orig.timepoint")
VlnPlot(cells, features = c("Id2","Id3"), pt.size = 0, assay = "RNA", layer = "counts")
VlnPlot(cells, features = "Prdm1", pt.size = 0, assay = "integrated")

VlnPlot(cells, features = c("Id2","Id3"), pt.size = 0, assay = "integrated")

