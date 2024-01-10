### trackplot prep script
require(tidyverse)
require(rtracklayer)
## source script from github
source("R/trackplot.R")
## TSS location per each gene (grep for specific gene)
ncbi_genes <- read.delim("data/annotations/mm39.ncbiRefSeq.gtf.gz", header = FALSE)
ncbi_genes <- separate(ncbi_genes, V9, into = c("gene_id", "transcript_id", "gene_name"),
                    sep = ";\\s*|;\\s*", extra = "drop", fill = "right") %>% 
                    #filter(gene_name == "exon_number 1") %>%
                    mutate(across(everything(), ~ gsub(".* ", "", .)))
## Test
ncbi_genes[ncbi_genes$gene_id == "Tbx21",]
colnames(ncbi_genes) <- c("chr", "library", "type", "start", "end", "score", "strand", "frame", "gene_id", "transcript_id", "gene_name")
## read GFF mm39
encode_gff <- read.delim("data/annotations/mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021.gff",
                         header = FALSE, comment.char = "#")
encode_gff <- separate(encode_gff, V9, into = c("id", "b_end", "b_start", "description", "feature_type"),
                       sep = ";", extra = "drop", fill = "right") %>%
                       mutate(across(everything(), ~ gsub(".*:", "", .)))
colnames(encode_gff) <- c("chr", "library", "type", "start", "end", "score", "strand", "frame", "id", "b_end", "b_start", "description", "feature_type")
##==========================## function wrapper ##==========================## 
track_gene <- function(
    gene = "Tbx21",
    reg_col = "#BEC446",
    ensembl_build = "mm39"
  ) {
  if (ensembl_build == "mm10") {
    bigWigs = c(
        "data/bigwig/GSE100738_NK.27-11b+.Sp.bw",
        "data/bigwig/GSE100738_NKT.Sp.bw",
        "data/bigwig/GSE100738_Tgd.Sp.bw",
        "data/bigwig/GSE100738_T.8.Nve.Sp.bw",
        "data/bigwig/GSE100738_T.4.Nve.Sp.bw"
        )
  } else if (ensembl_build == "mm39") {
    ## find gene transcript
    gene_results <- ncbi_genes[ncbi_genes$gene_id == gene & ncbi_genes$type == "transcript", ]
    gene_results <- gene_results[1,]
    ## bigWig files in mm39 (liftover)
    bigWigs = c(
        "data/mm39bw/GSE100738_NK.27-11b+.Sp_mm39.bw",
        "data/mm39bw/GSE100738_NKT.Sp_mm39.bw",
        "data/mm39bw/GSE100738_Tgd.Sp_mm39.bw",
        "data/mm39bw/GSE100738_T.8.Nve.Sp_mm39.bw",
        "data/mm39bw/GSE100738_T.4.Nve.Sp_mm39.bw"
        )
    ## define regions of interest
    reg_regions <- encode_gff[encode_gff$chr == sub("chr", "", gene_results$chr) & 
                              (encode_gff$start > as.numeric(gene_results$start)-10000) & 
                              (encode_gff$end < as.numeric(gene_results$end)+10000),]
  } else {
    message("error, wrong ensembl build selected")
  }
  # Extract the signal for your loci of interest ## T-bet on chr11:97,020,000-96,985,000
  track_data = track_extract(bigWigs = bigWigs,
                             loci = paste0(gene_results$chr,
                                           ":",
                                           if(gene_results$strand == "-") { 
                                             paste(gene_results$start)
                                           } else {
                                             paste(as.numeric(gene_results$start)-10000)},
                                           "-",
                                           if(gene_results$strand == "-") { 
                                             paste(as.numeric(gene_results$end)+10000)
                                           } else {
                                             paste(gene_results$end)}))# "chr11:97,000,897-97,020,157"
  ### add names to each dataset. 
  track_data = track_summarize(summary_list = track_data, 
                               condition = c("1 NK cells", "2 NKT cells", "3 gd T cells", "4 CD8 T cells", "5 CD4 T cells"), stat = "mean")
  ## Highlight regions of interest
  markregions <- reg_regions
  markregions$chr <- paste0("chr",markregions$chr)
  markregions %>% dplyr::select(chr, start, end, type) -> markregions
  colnames(markregions)[4] <- "name"
  ## JASPAR 2022
  cmd <- paste0(
    "bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/mm39/jaspar/JASPAR2022.bb -chrom=chr",reg_regions[reg_regions$type == "promoter","chr"], 
    " -start=", reg_regions[reg_regions$type == "promoter","start"],
    " -end=", reg_regions[reg_regions$type == "promoter","end"], 
    " output/data/", gene, ".bed"
  )

  ## run track plot
  track_plot(
    summary_list = track_data,
    draw_gene_track = TRUE,
    show_ideogram = TRUE,
    build = ensembl_build,
    regions = markregions,
    boxcol = reg_col,
    genename = gene_results$gene_symbol,
    gene_model = "data/annotations/refGene.gtf.gz",
    isGTF = TRUE
  )
}
