#!usr/bin/env Rscript
##
## load libraries
library(biomaRt)
library(Gviz)
library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(tidyverse)
library(TFBSTools)
## set databases
bm <- useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
## genes to querry
genes <- read.csv("data/gene_list.csv", header = FALSE)$V1
## bigWig ImmGEN files in mm39 (liftover)
## needs a lot of memory to hold on to all these files
bw <- c(
    "data/mm39bw/GSE100738_NK.27-11b+.Sp_mm39.bw",
    "data/mm39bw/GSE100738_NKT.Sp_mm39.bw",
    "data/mm39bw/GSE100738_Tgd.Sp_mm39.bw",
    "data/mm39bw/GSE100738_T.8.Nve.Sp_mm39.bw",
    "data/mm39bw/GSE100738_T.4.Nve.Sp_mm39.bw"
    )
##
## genes mm39
ncbi_genes <- read.delim("data/annotations/mm39.ncbiRefSeq.gtf.gz", header = FALSE)
ncbi_genes <- separate(ncbi_genes, V9, into = c("gene_id", "transcript_id", "gene_name"),
                    sep = ";\\s*|;\\s*", extra = "drop", fill = "right") %>% 
                    #filter(gene_name == "exon_number 1") %>%
                    mutate(across(everything(), ~ gsub(".* ", "", .)))
colnames(ncbi_genes) <- c("chr", "library", "type", "start", "end", "score", "strand", "frame", "gene_id", "transcript_id", "gene_name")

## regulatory regions mm39
encode_gff <- read.delim("data/annotations/mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021.gff",
                         header = FALSE, comment.char = "#")
encode_gff <- separate(encode_gff, V9, into = c("id", "b_end", "b_start", "description", "feature_type"),
                       sep = ";", extra = "drop", fill = "right") %>%
                       mutate(across(everything(), ~ gsub(".*:", "", .)))
colnames(encode_gff) <- c("chr", "library", "type", "start", "end", "score", "strand", "frame", "id", "b_end", "b_start", "description", "feature_type")
encode_gff$start <- as.numeric(encode_gff$start)
encode_gff$end <- as.numeric(encode_gff$end)
encode_gff$strand <- "*"
encode_gff$chr <- paste0("chr",encode_gff$chr)
encode_gff <- encode_gff[!grepl("\\.1$", encode_gff$chr), ]
encode_gff <- encode_gff[!grepl("\\.2$", encode_gff$chr), ]
## load bw files with rtracklayer
bw_file <- import.bw(bw[1], as = "GRanges") ## make this a loop or apply?
bw_fil2 <- import.bw(bw[2], as = "GRanges")
bw_fil3 <- import.bw(bw[3], as = "GRanges")
bw_fil4 <- import.bw(bw[4], as = "GRanges")
bw_fil5 <- import.bw(bw[5], as = "GRanges")

## get DNA string gene information from mm39
get_dna <- function(gene_name) {
    ## ncbi data 
    res <- ncbi_genes[ncbi_genes$gene_id == gene_name & ncbi_genes$type == "CDS",]
    res <- res[1,]
    ## get dna
    require(BSgenome.Mmusculus.UCSC.mm39) ## 600 Mb ## to get string of DNA)
    if (res$strand == "+") {
        dna <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm39, res$chr,
                              as.numeric(res$start) - 10000,
                              as.numeric(res$end)
                              )
    } else {
        dna <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm39, res$chr,
                              as.numeric(res$start),
                              as.numeric(res$end) + 10000
                              )
    }
    out <- list(res, dna)
    return(out)
}

#=======================# draw into Gviz tracks #=======================#
##
## use querry gene to generate list of sites with scores
draw_tracks <- function(q = genes, genome = "mm39") {
    if (genome == "mm39") {
        gen <- "mm39"
        message("starting run on mm39")
    } else {
        message("genome not supported, use mm39")
    }
    ## loop 
    for (gene in q) {
        message(paste("Beginning analysis of",gene))
        querry <- get_dna(gene)
        ## tracks parameters
        if (querry[[1]]$strand == "-" && gene != "Klrg1" && gene != "Klrb1c") {
            lims <- c(as.numeric(querry[[1]]$start) - 1000,
                      as.numeric(querry[[1]]$end) + 12000)
        } else if (querry[[1]]$strand == "-" && gene == "Klrb1c") {
           lims <- c(as.numeric(querry[[1]]$start) - 1000,
                      as.numeric(querry[[1]]$end) + 40000)
        } else if (querry[[1]]$strand == "-" && gene == "Klrg1") {
           lims <- c(as.numeric(querry[[1]]$start) - 1000,
                      as.numeric(querry[[1]]$end) + 20000)
        } else {
            lims <- c(as.numeric(querry[[1]]$start) - 12000,
                      as.numeric(querry[[1]]$end) + 1000)
        }
        print(querry)
        ## ATAC-seq data
        ## adjust ylim
        if (gene == "Nkg7") {
            ylims <- c(0, 5.5)
        } else {
            ylims <- c(0, 4.5)
        }
        ## make into loop of tracks?
        dtrack <- DataTrack(range = bw_file,
                    genome = gen,
                    type = "l",
                    chromosome = querry[[1]]$chr,
                    name = "NK",
                    start = lims[1], end = lims[2],
                    ylim = ylims
        )
        dtrac2 <- DataTrack(range = bw_fil2,
                    genome = gen,
                    type = "l",
                    chromosome = querry[[1]]$chr,
                    name = "NKT",
                    start = lims[1], end = lims[2],
                    ylim = ylims
        )
        dtrac3 <- DataTrack(range = bw_fil3,
                    genome = gen,
                    type = "l",
                    chromosome = querry[[1]]$chr,
                    name = "gdT",
                    start = lims[1], end = lims[2],
                    ylim = ylims
        )
        dtrac4 <- DataTrack(range = bw_fil4,
                    genome = gen,
                    type = "l",
                    chromosome = querry[[1]]$chr,
                    name = "CD8",
                    start = lims[1], end = lims[2],
                    ylim = ylims
        )
        dtrac5 <- DataTrack(range = bw_fil5,
                    genome = gen,
                    type = "l",
                    chromosome = querry[[1]]$chr,
                    name = "CD4",
                    start = lims[1], end = lims[2],
                    ylim = ylims
        )
        message("data tracks loaded")
        ## get mm39 from biomart for genecode
        biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    chromosome = querry[[1]]$chr,
                                    start = lims[1],
                                    end = lims[2],
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    transcriptAnnotation = "symbol"
                                    )
        ## other tracks
        gtrack <- GenomeAxisTrack()
        itrack <- IdeogramTrack(genome = "mm10", chromosome = querry[[1]]$chr)
        message("accesory tracks loaded")
        ## try
        try(
            ## CpG islands
            cpgIslands <- UcscTrack(genome = gen, 
                        chromosome = querry[[1]]$chr,
                        track = "cpgIslandExt",
                        table = "cpgIslandExt",
                        from = lims[1],
                        to = lims[2],
                        trackType = "AnnotationTrack",
                        start = "chromStart", end = "chromEnd",
                        id = "name", shape = "box", fill = "#006400",
                        name = "CpG Islands")
        )
        try(
            ## jaspar TF motifs
            jaspar_track <- UcscTrack(genome = gen,
                        chromosome = querry[[1]]$chr,
                        track = "JASPAR Transcription Factors",
                        table = "jaspar2022",
                        from = lims[1],
                        to = lims[2],
                        trackType = "AnnotationTrack",
                        start = "chromStart", end = "chromEnd",
                        score = "score",
                        id = "TFName", shape = "box", fill = "#006400",
                        name = "JASPAR",
                        featureAnnotation = "TFName")
        )
        ## filter jaspar
        filt_scores <- jaspar_track@dp@pars$score > 650 ### UCSC regular filter 400
        filt_range <- jaspar_track@range[filt_scores]
        try(
            filt_jaspar <- AnnotationTrack(range = filt_range,
                               name = "Jaspar 2022",
                               genome = gen, chromosome = querry[[1]]$chr,
                               id = filt_range$id,
                               shape = "box", fill = "#006400",
                               featureAnnotation = "id",
                               fontcolor.feature = 1)
        )
        try(
            ## regulatory build mm39
            reg_track <- AnnotationTrack(range = encode_gff, gen = gen,
                             chromosome = querry[[1]]$chr,
                             from = lims[1],
                             to = lims[2],
                             name = "Regulatory Build",
                             id = encode_gff$id,
                             #featureAnnotation = "id",
                             #fontcolor.feature = 1,
                             group = encode_gff$type ## groupAnnotation = "group" 
                             #feature = encode_gff$type
                             )
        )
        feature(reg_track) <- encode_gff$type
        if(!exists("cpgIslands")) {
            message("couldn't get CpG Islands")
            plist <- list(itrack,
                        gtrack,
                        biomTrack,
                        reg_track,
                        dtrack, dtrac2, dtrac3, dtrac4, dtrac5,
                        filt_jaspar)
        } else {
            plist <- list(itrack,
                        gtrack,
                        biomTrack,
                        reg_track,
                        dtrack, dtrac2, dtrac3, dtrac4, dtrac5,
                        filt_jaspar,
                        cpgIslands)
        }
        ## plot all tracks
        pdf(paste0("output/figures/tracks/gviz_",gene,".pdf"))
        plotTracks(plist,
                chromosome = querry[[1]]$chr,
                from = lims[1],
                to = lims[2]
           )
        dev.off()
        message(paste("Plot for", gene, "is ready!"))
    }
    message("all genes are done! Check results.")
}

## run function :D
draw_tracks()