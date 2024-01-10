#!usr/bin/env Rscript
##
#
# BiocManager::install("org.Mm.eg.db")
#
## load libraries
library(tidyverse)
library(org.Mm.eg.db)
library(KEGGREST)
## load data
beta_lvs <- read.table(file = "output/b_levels/results_mouse_beta_table_filt.csv", sep = ",", header = TRUE)
#
# Pentose-phosphate pathway ?
# KEGG
## Get KEGG pathways
pathways <- keggList("pathway", "mmu")
#
# Pull all genes for each pathway
pathway.codes <- gsub("path:", "", names(pathways))
genes.by.pathway <- sapply(pathway.codes,
                           try(
                             function(pwid){
                               pw <- keggGet(pwid)
                               if (is.null(pw[[1]]$GENE)) return(NA)
                               pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                               pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                               return(pw2)
                             }
                           )
)
head(genes.by.pathway)
##
genes.by.pathway <- genes.by.pathway[!is.na(genes.by.pathway)]


##
### run
data <- beta_lvs
genes <- data[data$pval < 0.05,]$gene
mat <- AnnotationDbi::select(org.Mm.eg.db, keys = genes, columns = c("ENTREZID","GENENAME"), keytype="SYMBOL")
## make gene list
geneList <- data[data$gene %in% mat$SYMBOL,]$pval
names(geneList) <- data[data$gene %in% mat$SYMBOL,]$gene
names(geneList) <- mat[match(names(geneList),mat$SYMBOL),]$ENTREZID
# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                               function(pathway) {
                                 pathway.genes <- genes.by.pathway[[pathway]]
                                 list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                                 list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                                 scores.in.pathway <- geneList[list.genes.in.pathway]
                                 scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                                 if (length(scores.in.pathway) > 0){
                                   p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
                                 } else{
                                   p.value <- NA
                                 }
                                 return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))
                               }
))
# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]

## plot libs
library(ggplot2)
outdat %>% filter(p.value < 0.05) %>% ggplot() +
  aes(x = reorder(pathway.name, -p.value), y = -log10(p.value)) + geom_col(fill = "#727d7e") + coord_flip() + theme_bw()
ggsave(filename = "output/figures/kegg_innate.pdf", width = 15, height = 8)
## write text
write.csv(x = outdat, file = paste0("output/data/kegg_","beta_filt.csv"))






