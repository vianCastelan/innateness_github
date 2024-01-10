#!/usr/bin/env Rscript

## load data
#beta_lvs <- read.table(file = "output/b_levels/results_mouse_beta_table.tsv", sep = "\t", header = TRUE)

## volcano_function
volcanoplot <- function(degtable, psig = 0.05, fcsig = 2, method = "MAST") {
  #' Make Volcanoplot
  #' 
  #' @description
  #'  `volcanoplot()` returns a `ggplot`  from values give in degtable generated from MAST or DESeq2.
  #' 
  #' # variables
  #' 
  #'  will require a `j` and a `genotype` variables to be define in the global environment.
  #' 
  #' # methods
  #' 
  #' @requirements 
  #' This functions requires ggplot2, ggrepel and maybe tidyverse.
  require(tidyverse)
  require(ggplot2)
  require(ggrepel)
  if(method == "MAST") {
    xval <- "avg_log2FC"
    yval <- "p_val_adj"
    labval <- "DE_status"
  } else if(method == "DESeq2") {
    xval <- "log2FoldChange"
    yval <- "padj"
    labval <- "DE_status"
  } else if(method == "beta") {
    xval <- "beta_norm"
    yval <- "pval"
    labval <- "DE_status"
  }
  ## def vline 
  vline = fcsig/2
  ### make DE_status column 
  dtb <- degtable
  dtb$genes <- if (method == "beta") dtb$gene else rownames(dtb)
  dtb$DE_status <- ifelse((dtb[,xval] >= fcsig) & (dtb[,yval] < psig), "UP",
                          ifelse((dtb[,xval] < fcsig*(-1)) & (dtb[,yval] < psig), "DOWN", "Not Significant"))
  ## ggplot 
  ggplot(dtb) + aes(x=get(xval), 
                    y=-log10(get(yval)), 
                    col=get(labval)) + geom_point(alpha = .4, size =2) + 
    scale_colour_manual(values=c("UP"="darkred","DOWN"="darkblue", "Not Significant"="gray")) +
    labs(title = paste("beta levels volcanoplot ","1000 MVG"), x = "beta levels", y = "log10 p-value") +
    geom_text_repel(data = dtb[dtb$DE_status != "Not Significant",],
                    aes(label = dtb[dtb$DE_status != "Not Significant",]$genes),
                    show.legend = FALSE, max.overlaps = 40, size = 8) +
    geom_vline(xintercept = c(-vline, 0, vline), color="darkgray") +
    theme_bw() + theme(legend.position = "none", text=element_text(size=18))
}

# make bi-exponential scale
biexp2_trans <- function(lim = 10000, decade.size = lim) {
    trans <- function(x) {
        ifelse(x <= lim,
            x,
            lim + decade.size * (suppressWarnings(log(x, 10)) -
                log(lim, 10))
        )
    }
    inv <- function(x) {
        ifelse(x <= lim,
            x,
            10^(((x - lim) / decade.size) + log(lim, 10))
        )
    }
    breaks <- function(x) {
        if (all(x <= lim)) {
            pretty_breaks()(x)
        } else if (all(x > lim)) {
            log_breaks(10)(x)
        } else {
            unique(c(
                pretty_breaks()(c(x[1], lim)),
                log_breaks(10)(c(lim, x[2]))
            ))
        }
    }
    trans_new(paste0("biexp-", format(lim)), trans, inv, breaks)
}

## volcano for RNA-seq
beta_lvs <- read.table(file = "output/b_levels/results_mouse_beta_table.csv", sep = ",", header = TRUE)
volcanoplot(beta_lvs[abs(beta_lvs$beta) > 0.2,], psig = 0.2,method = "beta") + scale_x_continuous(trans="biexp2") + theme(legend.position = "none")
ggsave("output/figures/beta_volcano.pdf", height = 7, width = 7)
rna <- beta_lvs
# ## load ATAC-seq data (filter levels close to 1)
beta_lvs <- read.table(file = "output/b_levels/results_mouse_beta_atac.tsv", sep = "\t", header = TRUE)
library(scales)
beta_lvs$beta_norm <- rescale(beta_lvs$beta, to=c(-2,2))
beta_lvs$beta_norm <- (beta_lvs$beta_norm - 1.26)*-1
beta_lvs <- beta_lvs[abs(beta_lvs$beta) > 0.2,]
volcanoplot(beta_lvs, psig = 0.05, fcsig = 0.2, method = "beta") + scale_x_continuous(trans="biexp2") + xlim(-2,2) +theme(legend.position = "none")
ggsave("output/figures/beta_volcano_atac.pdf", height = 7, width = 7)


## Dual VolcanoPlot 
comb <- merge(rna, beta_lvs, by = "gene")
head(comb)
colnames(comb)[5] <- "pval_rna"
colnames(comb)[10] <- "b_norm_rna"
colnames(comb)[14] <- "pval_atac"
colnames(comb)[16] <- "b_norm_atac"
comb$DE_status <- ifelse((comb[,"b_norm_rna"] >= 0.1) & (comb[,"pval_rna"] < 0.05), "UP",
                          ifelse((comb[,"b_norm_rna"] < 0.1*(-1)) & (comb[,"pval_rna"] < 0.05), "DOWN", "Not Significant"))
## plot 
comb %>% distinct(gene, .keep_all = TRUE) %>%
    ggplot() +
    aes(
        x = b_norm_rna,
        y = b_norm_atac,
        col = DE_status
    ) + geom_point(alpha = .4, size =2) + 
    scale_colour_manual(values=c("UP"="darkred","DOWN"="darkblue", "Not Significant"="gray")) +
    scale_x_continuous(trans="biexp2") +
    scale_y_continuous(trans="biexp2") + 
    geom_text_repel(data = comb %>% distinct(gene, .keep_all = TRUE) %>% filter(DE_status != "Not Significant"),
                    aes(label = (comb %>% distinct(gene, .keep_all = TRUE) %>% filter(DE_status != "Not Significant"))$gene),
                    show.legend = FALSE, max.overlaps = 20, size = 8) +
    xlim(-2,2) + theme_bw() + theme(legend.position = "none", text=element_text(size=18))
ggsave("output/figures/dual_rna_atac_volcanoplot.pdf", height = 7, width = 7)
