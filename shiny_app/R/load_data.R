# Bulk RNA-seq data
# --------------------------------------------------------------------------
count <- read.csv(file = "data/immgen_ULI_RNAseq.csv", header = TRUE, check.names = FALSE)
colnames(count)[1] <- "gene"

# Metadata
# --------------------------------------------------------------------------
metadata <- read.csv(file = "data/immgen_ULI_RNAseq_metadata.csv", header = TRUE)
rownames(metadata) <- metadata$names
metadata$X <- NULL

# Beta Table
# --------------------------------------------------------------------------
beta <- read.table(file = "data/results_mouse_beta_filt.tsv", sep = "\t", header = TRUE)