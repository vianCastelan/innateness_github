# boxplot_plot
plot_boxplot <- function(gene_to_plot, font_size = 32) {
    require(tidyverse)
    require(ggplot2)
    #require(patchwork)
    # labs 
    labx <- c("CD4","CD8","MAIT","gdT","iNKT","NK")
    # filter counts and col_data
    col_data <- metadata[metadata$tissue == "Spleen" & metadata$activation == "Naive", ]
    counts <- count[, c("gene", col_data$names)]
    # ranks for celltype
    col_data$rank <- c(rep(6,5), rep(5,3), rep(1,6),rep(2,4), rep(4,2), rep(3,2))
    merge(col_data,
        counts %>%
        filter(gene == gene_to_plot) %>%
        gather(key = "names", value = gene_to_plot)) %>%
    ggplot() +
    aes(x = reorder(all_conditions, rank),
      y = as.numeric(gene_to_plot), fill = all_conditions) +
    geom_boxplot() +
    theme_bw() + theme(text = element_text(size = font_size), 
                       legend.position = "none",
                       axis.text = element_text(size = font_size),
                       axis.text.x = element_text(angle=45, 
                                                  margin = margin(t = 25))
                       ) +
    ylab(substitute(italic(x)~"(CPM)", list(x=gene_to_plot))) + xlab(" ") + 
    ggtitle(substitute(italic(x)~"expression", list(x=gene_to_plot))) +
    scale_x_discrete(labels = labx)
}