library(optparse)
library(DropletUtils)

#https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html

option_list <- list(
  make_option(
    c("--raw_10x_dir"), type = "character", help = "Directory containing raw 10x counts"
  ),
  make_option(
    c("--filtered_10x_dir"), type = "character", help = "Directory containing filtered 10x counts"
  ),
  make_option(
    c("--filter_log"), default=NA, type = "character"
  ),
  make_option(
    c("--counts_nuclei"), type = "character", help = "Output file name for storing quality nuclei counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--counts_empty"), type = "character", help = "Output file name for storing empty droplets counts matrix (dgCMatrix format)."
  )
)


option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# Load counts for empty and all droplets
cells <- read10xCounts(opts$filtered_10x_dir, col.names=TRUE)
all_droplets <- read10xCounts(opts$raw_10x_dir, col.names=TRUE)

# Convert to matrix
matrix_cells <- counts(cells)
matrix_all <- counts(all_droplets)

# Empty droplets = all cells filtered out by 10x's cellranger
matrix_empty <- matrix_all[, setdiff(colnames(matrix_all), colnames(matrix_cells))]

# Print sizes
print(paste0("Empty matrix size: ", dim(matrix_empty)[2]))
print(paste0("Cell matrix size: ", dim(matrix_cells)[2]))

saveRDS(matrix_cells, file=opts$counts_nuclei)
saveRDS(matrix_empty, file=opts$counts_empty)

if (!is.na(opts$filter_log)) {
# Write # filtered to output file
line = paste("cellranger_cells", dim(matrix_all)[2], dim(matrix_cells)[2],  sep=",")
write(line,file=opts$filter_log,append=TRUE)

}

