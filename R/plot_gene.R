#!/usr/bin/env Rscript
### immgen dataset querry script
library(tidyverse)
library(yaml)
## read yaml 
config <- read_yaml("config.yaml") ## reads yaml as a list
if(version$major == 4 & version$minor > 0.3) {
  setwd(config$innate_dir$server)
  outdir <- config$innate_output$server
} else {
  setwd(config$innate_dir$local)
  outdir <- config$innate_output$local
}
cat("working directory is: ", getwd(), "\n")
## parse arguments
## Collect arguments
args <- commandArgs(TRUE)
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)
## Default setting when no all arguments passed or help needed
if("--help" %in% args | is.null(args$gene) | is.null(args$tissue)) {
  cat("
      The R Script arguments_section.R
      
      Mandatory arguments:
      --gene=text           - Use either MGI symbol for mouse genes or HGNC symbols for human genes.
      --tissue=text         - Either spleen, lung, liver or thymus.
      --help                - print this text

      Optional arguments:

      --mem=true/false      - Use to include CD4 and CD8 T effector populations.
      

      WARNING : This script requires gene and tissue input. 
      
      Example:
      Rscript R/plot_gene.R --gene=Cxcr6 --tissue=spleen \n
      or \n
      Rscript R/plot_gene.R --gene=Cxcr6 --tissue=spleen --mem=true
      \n")

  q(save="no")
}
cat("Gene selected : ", args$gene,"\n",sep="")
cat("Tissue selected : ", args$tissue,"\n",sep="")
if (args$mem == "true") cat("Memory populations will be included","\n", sep="")
## load data
immgen <- read.csv("data/countmatrix/GSE109125_filtered_genes.csv", header = TRUE, check.names = FALSE) ; colnames(immgen)[1] <- "gene"
col_data <- read.csv("data/metadata/immgen_ULI_RNAseq_metadata.csv", header = TRUE, row.names = 1)
## FILTER 
col_data %>% filter(
  tissue == "Spleen",
  #cell_type == "MAIT",
  activation == "Naive"
                    ) -> metadata
if (args$mem == "true") {
  metadata <- rbind(col_data %>% filter(names %in% c("T8.TE.LCMV.d7.Sp#1",
                                                     "T8.TE.LCMV.d7.Sp#2",
                                                     "T.4.Sp.aCD3+CD40.18hr#1",
                                                     "T.4.Sp.aCD3+CD40.18hr#2")),metadata)
  metadata$rank <- c(rep(4,2), rep(2,2), rep(8,5), rep(7,3), rep(1,6),rep(3,4), rep(6,2), rep(5,2))
} else {
  metadata$rank <- c(rep(6,5), rep(5,3), rep(1,6),rep(2,4), rep(4,2), rep(3,2))
}
## filter gene expression
immgen <- immgen[,c("gene",col_data$names)] 
## gene to plot
if (is.null(args$gene)) {
  gene_to_plot = "Nfatc1"
} else {
  gene_to_plot = args$gene
}
## plot 
merge(metadata, 
      immgen %>% 
        filter(gene == gene_to_plot) %>% 
        gather(key = "names", value = gene_to_plot)) %>% filter(all_conditions != "MAIT:Lung:Naive_Set") %>% 
  ggplot() + 
  aes(x = reorder(all_conditions, rank), 
      y = as.numeric(gene_to_plot), fill = all_conditions) + 
  geom_boxplot() + #geom_point() +
  theme_bw() + theme(text = element_text(size = 20), 
                     legend.position = "none", 
                     axis.text.x = element_text(angle = 90),
                     axis.title.y = element_text(face = "italic")
                     ) +
  ylab(paste(gene_to_plot, "(CPM)")) + xlab(" ") 
  #scale_fill_manual(values = c("red", "dark green", "purple", "light blue"))
ggsave(paste0("output/figures/",args$gene,".pdf"))
