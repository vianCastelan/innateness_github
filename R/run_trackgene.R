#!usr/bin/env Rscript

## run trackplot wrapper
source("R/trackplot_wrapper.R") ### this takes a while as SQL must be downloaded from several mirrors. Sometimes Ensembl site is unresponsive


# gene = "Klrb1c"
# pdf(paste0("tracks/",gene, "_tracks.pdf"))
# track_gene(gene = gene)
# dev.off()


genes <- c(
  "Klrg1",
  "Ccl5",
  "Nkg7",
  "Gzma",
  "Klrk1",
  "Il2rb",
  "Actb",
  "Prf1",
  "Klrb1c"
)

## make jaspar for these genes

## loop
for (gene in genes) {
  pdf(paste0("output/figures/tracks/",gene, "_tracks.pdf"))
  track_gene(gene = gene, 
    ensembl_build = "mm39"
  )
  dev.off()
}
