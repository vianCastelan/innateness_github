#!/usr/bin/env Rscript

## plot top 20 b-levels for innateness.

library(tidyverse)
library(ggplot2)
library(cowplot)
## parameters 
colors <- read.csv("colors.csv")


### Beta levels plots 
beta_lvs <- read.table(file = "output/b_levels/results_mouse_beta_table_filt.csv", sep = ",", header = TRUE)
human <- readxl::read_xls(path = "data/human/41467_2019_8604_MOESM6_ESM.xls", range = 'B3:I19934', col_names = TRUE)
colnames(human)[3] <- "beta"

## make 2 plots 
top_mouse <- rbind(
  ## innate mouse
  beta_lvs %>%
  top_n(n = 20, beta_norm),
  ## adaptive mouse
  beta_lvs %>%
  top_n(n = 20, -beta_norm)
) %>% ggplot() +
  aes(x = reorder(gene, beta_norm),
      y = beta_norm,
      fill = c(rep(colors$color[8], 20), rep(colors$color[7], 20))) +
  geom_col(col = "black") +
  coord_flip() + theme_classic() + theme(legend.position = "none")

top_human <- rbind(
  ## innate human
  human %>%
  top_n(n = 20, beta),
  ## adaptive human
  human %>%
  top_n(n = 20, -beta)
) %>% ggplot() +
  aes(x = reorder(GENE_NAME, beta),
      y = beta,
      fill = c(rep(colors$color[8], 20), rep(colors$color[7], 20))) +
  geom_col(col = "black") +
  coord_flip() + theme_classic() + theme(legend.position = "none")
## plot grid
plot_grid(top_mouse, top_human)
ggsave("output/figures/beta_human_mouse.pdf", width = 8, height = 8)
