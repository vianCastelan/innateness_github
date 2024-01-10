#!/usr/bin/env Rscript
# libraries
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
cellcol <- read.csv("colors.csv")
#Reading beta table output from linear mixed models
beta_th1 <- read.table(file = "output/b_levels/Th1_beta_table.tsv", header = TRUE, sep = '\t')
beta_th2 <- read.table(file = "output/b_levels/Th2_beta_table.tsv", header = TRUE, sep = '\t')
beta_th17 <- read.table(file = "output/b_levels/Th17_beta_table.tsv", header = TRUE, sep = '\t')


## filter
beta_th1 <- na.omit(beta_th1[beta_th1$pval < 0.01 & abs(beta_th1$beta) > 50,])
beta_th2 <- na.omit(beta_th2[beta_th2$pval < 0.01 & abs(beta_th2$beta) > 50,])
beta_th17 <- na.omit(beta_th17[beta_th17$pval < 0.01 & abs(beta_th17$beta) > 10,])

## VennDiagram
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(
    beta_th1$gene,
    beta_th2$gene,
    beta_th17$gene
  ),
  category.names = c("Th1" , "Th2" , "Th17"),
  filename = NULL, #'output/figures/Th_innate_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
) -> vnobj
ggsave(vnobj, file = "output/figures/Th_innate_venn.pdf")

## add 1000 MVG models 
beta <- read.table("output/b_levels/results_mouse_beta_table_filt.csv", header = TRUE, sep = ",")
## new Venn
myCol <- brewer.pal(4, "Pastel2")
venn.diagram(
  x = list(
    beta_th1$gene,
    beta_th2$gene,
    beta_th17$gene,
    beta$gene
  ),
  category.names = c("Th1" , "Th2" , "Th17", "Innate"),
  filename = NULL, #'output/figures/Th_vs_innate_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  #lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #rotation = 1
) -> vnobj
ggsave(vnobj, file = "output/figures/Th_vs_innate_venn.pdf")

## Reading innateness scores applied to data
innateness_table_score1 <- read.table("output/b_scores/Th1_innateness.tsv", header = TRUE, sep = "\t")
innateness_table_score2 <- read.table("output/b_scores/Th2_innateness.tsv", header = TRUE, sep = "\t")
innateness_table_score17 <- read.table("output/b_scores/Th17_innateness.tsv", header = TRUE, sep = "\t")

## plot against each other
innate <- innateness_table_score1[,-1]
innate$th1 <- innateness_table_score1$innateness_score
innate$th2 <- innateness_table_score2$innateness_score
innate$th17 <- innateness_table_score17$innateness_score
innate$tissue <- ifelse(
  grepl("Sp",innate$cell_type),"Spleen", ifelse(
    grepl("Lv",innate$cell_type), "Liver", ifelse(
      grepl("Lu",innate$cell_type), "Lung", "Small Intestine"
    )
  )
)

## compare function
compare <- function(df,name1,name2) {
    x = df[,name1]
    y = df[,name2]
    my_plot <- ggplot(df, aes(x, y, color=as.factor(cell_type),shape = factor(tissue))) +
    geom_point(size=3) +
    scale_shape_manual(name = "Tissue", values = c(17,15,18,16))+
    ggtitle(paste(name1,"vs",name2,sep = " ")) +
    labs(color = "Cell type") + theme_bw() +
    xlab(paste(name1,"innateness score",sep = " ")) +
    ylab(paste(name2,"innateness score",sep = " "))
  return(my_plot)
}

p1 <- compare(innate,"th1","th17")
p2 <- compare(innate,"th1", "th2")
p3 <- compare(innate,"th2","th17")

plot_grid(p1,p2,p3, nrow =1)
ggsave(filename = "output/figures/th_applied.pdf", width = 16, height = 4)
# #Function to plot volcano plots
# volcano <- function(df,name){
#   my_plot<-ggplot(df,aes(x=(beta),y=-log10(pval),text=gene)) + geom_point() +
#     ggtitle(paste(name,"genes",sep=" ")) +
#     xlab(paste0("Innateness score")) +
#     ylab(paste0("-log10(pval)"))
#   top_peaks <- df[with(df, order(beta, pval)),][1:15,]
#   top_peaks <- rbind(top_peaks, df[with(df, order(-beta, pval)),][1:15,])
#   my_plot <- my_plot + geom_label_repel(data = top_peaks,aes(x=(beta),y=-log10(pval),label=gene),
#                               box.padding   = 0.08,
#                               point.padding = 0.3,
#                               segment.color = 'grey50') + theme_classic()
#   return(my_plot)
# }

# volcano(beta_129,"129 PCA")

# #Function to plot box plots
# box <- function(innateness_t,name){
#   my_plot <- ggplot(innateness_t, aes(x = reorder(cell_type,innateness_score), y = innateness_score, fill=as.factor(cell_type))) +
#     #scale_y_reverse()+
#     geom_boxplot() +
#     ggtitle(paste(name,"innatenes",sep = " ")) +
#     xlab(paste0("Cell type")) +
#     ylab(paste0("Innateness score"))+
#     coord_flip()+
#     theme_bw()+
#     theme(legend.position = "none")
#   return(my_plot)
# }

# my_plot1000<-box(beta_1000,"1000")
# my_plot129<-box(beta_129,"129")

# #Function to plot violin plots
# violin <- function(innatenes_t, name){
#   my_plot<-ggplot(innatenes_t, aes(x = reorder(cell_type,innateness_score), y = innateness_score, fill=as.factor(cell_type)))+
#     #scale_y_reverse()+
#     geom_violin() +
#     ggtitle(paste(name,"innatenes",sep = " ")) +
#     xlab(paste0("Cell type")) +
#     ylab(paste0("Innateness score"))+
#     coord_flip()+
#     theme_bw()+
#     theme(legend.position = "none")
#   return(my_plot)
# }

# #Functions to tabule bottom and top genes
# table_head <- function(df){
#   my_head <- df%>%select(gene,beta)%>%arrange(beta)%>%tail(20)
#   top <- tableGrob(my_head,
#                    theme=ttheme_default(base_size = 7))
#   return(top)
# }

# table_tail <- function(df){
#   my_tail <- df%>%select(gene,beta)%>%arrange(beta)%>%head(20)
#   bottom <- tableGrob(my_tail,
#                       theme=ttheme_default(base_size = 7))
#   return(bottom)
# }

# compare <- function(df,x,y,name1,name2) {
#     my_plot <- ggplot(df, aes(x, y, color=as.factor(cell_type),shape = factor(Tissue))) +
#     geom_point(size=3) +
#     scale_shape_manual(name = "Tissue", values = c(17,15,18,16))+
#     ggtitle(paste(name1,"vs",name2,sep = " ")) +
#     labs(color = "Cell type") + theme_bw() +
#     xlab(paste(name1,"innateness score",sep = " ")) +
#     ylab(paste(name2,"innateness score",sep = " "))
#   return(my_plot)
# }

# beta_tables <- list(beta_table1, beta_table2, beta_table17)
# innateness_tables <- list(innateness_table_score1, innateness_table_score2, innateness_table_score17)
# names <- list("Th1","Th2","Th17")
# volcano_plots <- list()
# box_plots <- list()
# violin_plots <- list()
# tops <- list()
# bottoms <- list()

# for(i in 1:3){
#   volcano_plots[[i]] <- volcano(beta_tables[[i]],names[[i]])
#   box_plots[[i]] <- box(innateness_tables[[i]],names[[i]])
#   violin_plots[[i]] <- violin(innateness_tables[[i]],names[[i]])
#   tops[[i]] <- table_head(beta_tables[[i]])
#   bottoms[[i]] <- table_tail(beta_tables[[i]])
# }

# pdf("plot_comparison_Th1_Th2_Th17_new_models.pdf")
# grid.arrange(volcano_plots[[1]], volcano_plots[[2]], volcano_plots[[3]], nrow = 1)
# grid.arrange(box_plots[[1]], box_plots[[2]], box_plots[[3]], nrow = 1)
# grid.arrange(violin_plots[[1]], violin_plots[[2]], violin_plots[[3]], nrow = 1)
# grid.arrange(tops[[1]], tops[[2]], tops[[3]], nrow = 1, 
#              top = textGrob("Th1, Th2, Th17 top innateness score genes",
#              gp = gpar(fontsize=20,font=3), x = 0.5, hjust = 0.5))
# grid.arrange(bottoms[[1]], bottoms[[2]], bottoms[[3]], nrow=1,
#              top = textGrob("Th1, Th2, Th17 bottom innateness score genes",
#              gp = gpar(fontsize=20,font=3), x = 0.5, hjust = 0.5))
# compare(all,all$Th1,all$Th17,"Th1","Th17")
# compare(all,all$Th1,all$Th2,"Th1","Th2")
# compare(all,all$Th2,all$Th17,"Th2","Th17")
# dev.off()