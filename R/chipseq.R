#!usr/bin/env Rscript
##
library(GEOquery)
## geo libraries
# GSE101593 CD4+ T cells SMAD4
# GSE55834 H3K4me3, H3K27me3, H3K4me3 NK cells (this is on mm9)
# GSE50129 NK, CD8+ T cell RUNX3
# GSE164906 multiple cells: CD4+, CD8+ NK cells ATAC and RNA... this is not very useful...
# GSE151637 early effector CD8+ T cells RNAseq
# GSE135533 CD8+ T cells SMAD4
# GSE58775 CD4+ T cell H3K4me3 and H3K27me3
# GSE75724 CD4+ T cells Blimp1
glist <- c("GSE101593", ## CD4+ T cell SMAD4, this one is no good... 
           "GSE55834", ## NK histones
           "GSE50129", ## NK Runx3
           "GSE135533", ## CD8+ T cell SMAD4
           "GSE58775", ## CD4+ H3K4me3 and H3K27me3, this failes to download locally 
           "GSE75724" ## CD4+ T cells Blimp1
           )

for (geo in glist) {
      filePaths <- getGEOSuppFiles(geo, baseDir = "data")
      tryCatch(untar(tarfile =  rownames(filePaths),
            exdir = paste0("data/", geo)),
            error = function(e) {
                  message(paste("An error occured when trying to un-tar",geo))
            })
      cat(geo, " done", "\n")
}

