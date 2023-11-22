library(Seurat)
library(DropletUtils)

# Parse args
args <- commandArgs(trailingOnly=TRUE)

# Read the object
seurat_object <- readRDS(args[1])

# Write to folder in 10X format
write10xCounts(args[2], seurat_object@assays$RNA@counts, version='3')
