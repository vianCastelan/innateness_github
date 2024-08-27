# IPA functions and plot functions
# Gabriel Ascui
# Version 10NOV20 || updated 18APR21
# requiered libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
########################################################  Parameters ########################################################
parameters <- list(
  zscore_sig = 0.5,
  pscore_sig = 0.05,
  lab_z_sig = 2,
  lab_p_sig = 0.01,
  RU = "RU",
  LU = "LU",
  colors = data.frame(
    colour = c(
    "#F68282",
    "#31C53F",
    "#1FA195",
    "#B95FBB",
    "#D4D915",
    "#28CECA",
    "#ff9a36",
    "#2FF18B",
    "#aeadb3",
    "#faf4cf",
    "#CCB1F1",
    "#25aff5",
    "#A4DFF2",
    "#4B4BF7",
    "#AC8F14",
    "#E6C122"
  ))
)
# cluster names
cnms <- c("1","2")
# Functions
################################            IPA analysis excel reader               #########################################
# canonical <- read_xls("27jan20/13PCs/DE_MAST/IPA_mapped/canpath_mait_10vs8.xls", skip = 1)
#################################  Pre-process IPA xls output for downstream plotting  ######################################
pp_canonical <- function(can = canonical, z_sig = 1, p_sig = 0.05){
  colnames(can) <- c("IPA", "log2p_val", "Ratio", "z_score", "molecules")
  can$n_genes <- sapply(can[,5], function(x) str_count(x, ",") +1)
  can$z_score[is.na(can$z_score)] <- 0
  # significance (updated 10NOV20)
  can$Significance <- "NS"
  can$Significance[(can$z_score > z_sig & can$log2p_val > -log2(p_sig))] <- "RU" # right-up
  can$Significance[(can$z_score < -(z_sig) & can$log2p_val > -log2(p_sig))] <- "LU" # left-up
  can$Significance <- factor(can$Significance, levels=c("NS", "RU", "LU")) 
  # next update, add optional z_score == 0 filtering
  return(can)
}
# Plotting functions
################################ bubble_plot function for IPA analysis presentation #########################################
bubble_plot <- function(can, bubble_size = 12){
  can %>% ggplot(aes(x = Ratio,
                         y = z_score,
                         color = log2p_val,
                         size = n_genes)) +
  geom_point(col="grey", alpha=.6) +
  geom_point(data = can %>% filter(abs(z_score) > parameters$lab_z_sig, Ratio > 0.2)) +
  geom_text(data = can %>% filter(abs(z_score) > parameters$lab_z_sig, log2p_val > -log2(parameters$lab_p_sig)), aes(label = IPA),
            nudge_y = 0.2, nudge_x = 0.01,
            check_overlap = TRUE, size = 3, color = "black"
  ) + scale_color_viridis_c() + theme_bw() + scale_size(range = c(.1, bubble_size), name="n_genes")
}
################################  can_plot function for IPA analysis presentation   #########################################
can_plot <- function(can = canonical, lb = path_labels, lgnd = FALSE, max_z = 3, max_p = 15) {
  if (!exists("path_labels")) {
    path_labels = subset(can, log2p_val > -log2(parameters$lab_p_sig) & abs(z_score) > parameters$lab_z_sig)
  }
  p = ggplot(can,aes(x = z_score,
                     y = log2p_val,
                     fill = Significance)) +
    xlim(-max_z,max_z) + ## set max_z
    ylim(0, max_p) +     ## set max_p
    scale_fill_manual(name = "Pathway",
                      labels = c("NS", parameters$RU,parameters$LU),
                      values=c(NS="grey", RU=parameters$colors[cnms[2],], LU=parameters$colors[cnms[1],])) +
    geom_hline(yintercept = -log2(parameters$pscore_sig), linetype="dashed", color = "grey") +
    geom_vline(xintercept = c(-(parameters$zscore_sig),parameters$zscore_sig), linetype = "dashed", color = "grey") +
    geom_point(shape = 21, size = 4) +
    geom_text_repel(data = path_labels, ## next update, make these values as variables
                    aes(label=IPA),
                    col = "black",
                    size = 4,
                    force =8,
                    box.padding = 0.5) +
    theme_bw(base_size = 14)
  if(lgnd == FALSE) {
    return(p + theme(legend.position = "none"))
  } else {
    return(p)
  }
}
################################  barcan_plot function for IPA analysis    #########################################
barcan_plot <- function(can = canonical, cluster = 1) {
  p =  ggplot(can, aes(x=reorder(IPA, log2p_val),
                       y=log2p_val)) +
    geom_col(fill=parameters$colors[cnms[cluster],],col="black") + coord_flip() + theme_classic()
  return(p)
}


################################              Other Functions              #########################################
# Function to convert strings from ALL CAPS to First Cap
to_first_cap <- function(s) {
  s <- tolower(s)
  s <- strsplit(s, " ")[[1]]
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2), sep="", collapse=" ")
  return(s)
}

# Function to remove content between parentheses
remove_parentheses <- function(s) {
  return(gsub("\\(.*?\\)", "", s))
}

# Function to process the input string
process_string <- function(input) {
  # Split the input string by comma
  parts <- strsplit(input, ",")[[1]]
  # Apply the transformations
  parts <- sapply(parts, function(x) to_first_cap(remove_parentheses(trimws(x))))
  parts <- unname(parts)
  return(parts)
}


