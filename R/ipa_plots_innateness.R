############ Script for IPA analysis on Innateness Mouse 
### Innateness Levels by Vianka Cedillo Castelan
## Gabriel Ascui 
## 18APR21
#
#
####==================#### source functions ####==================#### 
source("R/IPA_functions.R")  ## tidyverse, ggplot and ggrepel sourced here
# parameters 
parameters <- list(
  zscore_sig = 0.5,
  lab_z_sig = 2,
  lab_p_sig = 0.01,
  pscore_sig = 0.05,
  RU = "Innateness",
  LU = "Adaptiveness",
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
    )
  )
)
# group color select 
cnms <- c("1", "3") ### select color according to Rcolor_brewer
####==================####  load IPA data   ####==================#### 
canonical <- read.delim(file = "output/ipa/canonical_pathways_all.txt", skip = 2, check.names = FALSE)
canonical <- canonical[,-6]
canonical <- pp_canonical(canonical)
## Volcano Plot
#path_labels <- subset(canonical, log2p_val > 14 & abs(z_score) > 2)  ## this is not necessary to define anymore
canonical %>% filter(z_score != 0) %>% can_plot(lgnd = FALSE) + ylim(0,25)
ggsave(last_plot(), filename = "output/figures/ipa_innateness_adaptiveness.pdf", device = "pdf", width = 7, height = 6, dpi = 300)
####==================####    IPA plots Innateness Up    ####==================#### 
## load and pre-processing
canonical <- read.delim(file = "output/ipa/pos/canonical_pathways.txt", skip = 2, check.names = FALSE)
canonical <- canonical[,-6]
canonical <- pp_canonical(canonical)
path_labels <- subset(canonical, log2p_val > 8  & abs(z_score) > 1) ## specific path_labels 
## plot 
canonical %>% filter(z_score != 0) %>% can_plot(lgnd = FALSE) + ylim(0,25) + xlim(-6,6)
ggsave(last_plot(), filename = "output/figures/ipa_innateness_up_volcano.pdf", device = "pdf", width = 6, height = 6, dpi = 300)
canonical %>% filter(z_score != 0) %>% bubble_plot()
ggsave(last_plot(), filename = "output/figures/ipa_innateness_up_bubble.pdf", device = "pdf", width = 6, height = 6, dpi = 300)
canonical %>% filter(z_score != 0,
                     log2p_val > -log2(0.01), 
                     z_score > 0,
                     !grepl("B Cell", IPA)) %>% barcan_plot()
ggsave(last_plot(), filename = "output/figures/ipa_innateness_up_barcan.pdf", device = "pdf", width = 10, height = 6, dpi = 300)
## upstream regulators 
upstream <- read.delim(file = "output/ipa/pos/upsteam_regulators.txt", skip = 2, check.names = FALSE)
# pp_upstream
colnames(upstream) <- c("regulator", "Exp_FC", "type", "predicted_state", "z_score", "flags", "p_val", "targets", "mechanism")
# significance (updated 10NOV20)
upstream$Significance <- "NS"
upstream$Significance[(upstream$z_score > parameters$zscore_sig & upstream$p_val < parameters$lab_p_sig)] <- "RU" # right-up
upstream$Significance[(upstream$z_score < -(parameters$zscore_sig) & upstream$p_val < parameters$lab_p_sig)] <- "LU" # left-up
upstream$Significance <- factor(upstream$Significance, levels=c("NS", "RU", "LU")) 
## upstream_plot 
upstream %>% filter(type == "transcription regulator") %>% 
  ggplot() +
  aes(x = z_score, y = -log10(p_val), fill = Significance) + 
  geom_vline(xintercept = c(-(parameters$zscore_sig),parameters$zscore_sig), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(parameters$pscore_sig), linetype="dashed", color = "grey") +
  geom_point(shape = 21, size = 4) + 
  scale_fill_manual(name = "Programing",
                    labels = c("NS", parameters$RU,parameters$LU),
                    values=c(NS="grey", RU=parameters$colors[cnms[2],], LU=parameters$colors[cnms[1],])) +
  geom_text_repel(data = subset(upstream, p_val < parameters$lab_p_sig & abs(z_score) > parameters$zscore_sig & type == "transcription regulator"), ## next update, make these values as variables
                  aes(label=regulator),
                  col = "black",
                  size = 4,
                  force =8,
                  box.padding = 0.5) +
  theme_bw()
ggsave(last_plot(), filename = "output/figures/ipa_innateness_upstream_TFs.pdf", device = "pdf", width = 6, height = 5, dpi = 300)
####==================####    IPA plots Innateness Up    ####==================#### 
## load and pre-processing
#canonical <- read.delim(file = "output/ipa/neg/canonical_pathways.txt", skip = 2, check.names = FALSE)
#canonical <- canonical[,-6]
#canonical <- pp_canonical(canonical)
#path_labels <- subset(canonical, log2p_val > 8  & abs(z_score) > 1) ## specific path_labels 
## plot 
canonical %>% filter(z_score != 0) %>% can_plot(lgnd = FALSE) + ylim(0,30) + xlim(-6,6)
ggsave(last_plot(), filename = "output/figures/ipa_adaptiveness_up_volcano.pdf", device = "pdf", width = 6, height = 6, dpi = 300)
canonical %>% filter(z_score != 0) %>% bubble_plot()
ggsave(last_plot(), filename = "output/figures/ipa_adaptiveness_up_bubble.pdf", device = "pdf", width = 6, height = 6, dpi = 300)
canonical %>% filter(z_score != 0,
                     log2p_val > -log2(0.01), 
                     z_score > 0,
                     !grepl("B Cell", IPA)) %>% barcan_plot(cluster = 1)
ggsave(last_plot(), filename = "output/figures/ipa_adaptiveness_up_barcan.pdf", device = "pdf", width = 10, height = 6, dpi = 300)
####==================####    IPA plots Adaptiveness Down   ####==================#### 
canonical %>% filter(z_score != 0,
                     log2p_val > -log2(0.01), 
                     z_score < 0,
                     !grepl("B Cell", IPA)) %>% barcan_plot(cluster = 2)
ggsave(last_plot(), filename = "output/figures/ipa_adaptiveness_up_barcan.pdf", device = "pdf", width = 10, height = 6, dpi = 300)
