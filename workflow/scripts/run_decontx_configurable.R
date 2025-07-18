options(stringsAsFactors=FALSE)

library(optparse)
library(celda)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(yaml)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in raw counts matrices (dgCMatrix) for nuclei and empty droplets
# and runs decontx to generate decontaminated counts
# This version is configurable via YAML config files for different tissue types
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


option_list <- list(
  make_option(
    c("--counts_nuclei"), type = "character", help = "Input file name for storing quality nuclei counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--counts_empty"), type = "character", help = "Input file name for storing empty droplets counts matrix (dgCMatrix format)."
  ),
  make_option(
    c("--estimate_delta"), action='store_true', help = "Option to let decontX try to estimate delta (default FALSE)."
  ),
  make_option(
    c("--delta_first"), type = "numeric", help = "First delta value for decontX.\nCan be used whether --estimate_delta is TRUE or FALSE (see docs).", default=10
  ),
  make_option(
    c("--delta_second"), type = "numeric", help = "Second delta value for decontX.\nCan be used whether --estimate_delta is TRUE or FALSE (see docs).", default=30
  ),
  make_option(
    c("--clusters"), type = "character", help = "File mapping barcodes to clusters."
  ),
  make_option(
    c("--max_contamination"), type = "numeric", help = "Maximum contamination fraction for barcodes to return."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for decontX results."
  ),
  make_option(
    c("--step_name"), type = "character", help = "Name for step in filter log (optional)", default = NA
  ),
  make_option(
    c("--filter_log"), type = "character", help = "Filter log (optional)", default = NA
  ),
  make_option(
    c("--config"), type = "character", help = "Path to YAML configuration file for tissue-specific parameters"
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# Load configuration file
if (!is.null(opts$config) && file.exists(opts$config)) {
  config <- yaml.load_file(opts$config)
  cat("Loaded configuration from:", opts$config, "\n")
} else {
  stop("Configuration file not found or not specified: ", opts$config)
}

# Load marker genes configuration if specified
marker_config <- NULL
if (!is.null(config$marker_genes_config) && file.exists(config$marker_genes_config)) {
  marker_config <- yaml.load_file(config$marker_genes_config)
  cat("Loaded marker genes configuration from:", config$marker_genes_config, "\n")
}

# ### Testing
# opts = list()
# opts$counts_nuclei = "results/multiome/decontx/Sample_5124-NM-1-hg38/counts_nuclei.rds"
# opts$counts_empty = "results/multiome/decontx/Sample_5124-NM-1-hg38/counts_empty.rds"
# opts$clusters = "results/multiome/seurat_prelim/Sample_5124-NM-1-hg38/seurat_clusters.csv"
# opts$outdir = "results/multiome/decontx/Sample_5124-NM-1-hg38"
# opts$max_contamination = 0.3


### Read data
counts_nuclei = readRDS(opts$counts_nuclei)
counts_empty = readRDS(opts$counts_empty)
clusters = read.csv(opts$clusters)


# Filter counts_nuclei to only those found in clusters
counts_nuclei = counts_nuclei[, colnames(counts_nuclei) %in% clusters$barcode_gex]

### Message
cat("START: Training decontX model using", ncol(counts_nuclei), "nuclei and", ncol(counts_empty), "empty droplets.\n")

### Convert to SingleCellExperiment object
cat("Converting to SingleCellExperiment object.\n")
x_nuclei = SingleCellExperiment(assays = list(counts = counts_nuclei))
x_empty = SingleCellExperiment(assays = list(counts = counts_empty))

### Make sure cluster labels are in right order
row.names(clusters) = clusters$barcode_gex
z = clusters[colnames(counts_nuclei),"clusters"]

### Run decontx -- unclear if the background option is actually being used here
# Bug was supposedly fixed in dev branch: https://github.com/campbio/celda/issues/355
cat("Running decontX with the following parameters:\n")
cat("Estimate delta: ", opts$estimate_delta, "\n")
cat("Delta: c(", opts$delta_first, ", ", opts$delta_second, ")\n") 

results <- decontX(x = x_nuclei,
  z = z,
  background = x_empty,
  delta=c(opts$delta_first, opts$delta_second),
  estimateDelta = opts$estimate_delta,
  seed = 8675)


### Plot distribution of decontx values
contamination_estimates <- results$decontX_contamination
names(contamination_estimates) <- colnames(results)
contamination_estimates <- data.frame(contamination_estimates)

cells_removed_string <- paste0("Total cells: ", length(results$decontX_contamination), ", discarded: ", sum(results$decontX_contamination > opts$max_contamination))

print(head(contamination_estimates))
p <- ggplot(contamination_estimates) + 
        geom_histogram(aes(x=contamination_estimates)) +
        geom_vline(xintercept = opts$max_contamination, color='red') + 
        ggtitle(cells_removed_string)

ggsave(file.path(opts$outdir, "contamination_distribution.png"), height=4, width=6)

p <- ggplot(contamination_estimates) + 
        stat_ecdf(aes(x=contamination_estimates)) +
        geom_vline(xintercept = opts$max_contamination, color='red') + 
        ggtitle(cells_removed_string)

ggsave(file.path(opts$outdir, "contamination_ecdf.png"), height=4, width=6)


p <- ggplot(contamination_estimates) + 
        geom_density(aes(x=contamination_estimates)) +
        geom_vline(xintercept = opts$max_contamination, color='red') + 
        ggtitle(cells_removed_string)

ggsave(file.path(opts$outdir, "contamination_density.png"), height=4, width=6)


### UMAP plots
cat("Generating UMAP plots.\n")
umap <- reducedDim(results, "decontX_UMAP")

png(file.path(opts$outdir, "umap_clusters.png"))
plotDimReduceCluster(x = results$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])
dev.off()

png(file.path(opts$outdir, "umap_contamination.png"))
plotDecontXContamination(results)
dev.off()

# Get marker genes based on tissue type configuration
tissue_type <- config$tissue_type %||% "general"
known_markers <- c()

# Add general markers if available
if (!is.null(marker_config$known_markers)) {
  known_markers <- c(known_markers, marker_config$known_markers)
}

# Add tissue-specific markers if available
if (!is.null(marker_config$cell_type_markers[[tissue_type]])) {
  tissue_markers <- unlist(marker_config$cell_type_markers[[tissue_type]], use.names = FALSE)
  known_markers <- c(known_markers, tissue_markers)
} else if (!is.null(marker_config$general)) {
  # Fallback to general markers if tissue-specific ones aren't available
  general_markers <- unlist(marker_config$general, use.names = FALSE)
  known_markers <- c(known_markers, general_markers)
}

# Remove duplicates and ensure we have some markers
known_markers <- unique(known_markers)

# If no markers found in config, fall back to hardcoded pancreatic islet markers
if (length(known_markers) == 0) {
  cat("No markers found in configuration, using default pancreatic islet markers.\n")
  known_markers = c(
        "ENSG00000115263", # GCG
        "ENSG00000254647", # INS
        "ENSG00000171345", # KRT19
        "ENSG00000108849", # PPY
        "ENSG00000204983", # PRSS1
        "ENSG00000081237", # PTPRC
        "ENSG00000143248", # RGS5
        "ENSG00000232995", # RGS5
        "ENSG00000135094", # SDS
        "ENSG00000157005", # SST
        "ENSG00000110799"  # VWF
  )
} else {
  cat("Using", length(known_markers), "marker genes from configuration for tissue type:", tissue_type, "\n")
}

results <- logNormCounts(results)
png(file.path(opts$outdir, "umap_marker_genes_pre.png"))
plotDimReduceFeature(as.matrix(logcounts(results)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = known_markers,
    exactMatch = TRUE,
    useAssay = "counts"
  )
dev.off()

png(file.path(opts$outdir, "umap_marker_genes_post.png"))
plotDimReduceFeature(as.matrix(logcounts(results)),
    dim1 = umap[, 1],
    dim2 = umap[, 2],
    features = known_markers,
    exactMatch = TRUE,
    useAssay = "decontXcounts"
  )
dev.off()

# Create tissue-specific marker lists for barplots
markers <- list()

if (!is.null(marker_config$cell_type_markers[[tissue_type]])) {
  # Use configured markers for the specific tissue type
  tissue_marker_config <- marker_config$cell_type_markers[[tissue_type]]
  for (cell_type in names(tissue_marker_config)) {
    markers[[cell_type]] <- tissue_marker_config[[cell_type]]
  }
} else {
  # Fallback to default pancreatic islet markers
  cat("Using default pancreatic islet marker configuration for barplots.\n")
  markers <- list(
    beta_INS = c("ENSG00000254647"),
    alpha_GCG = c("ENSG00000115263"),
    delta_SST = c("ENSG00000157005"),
    gamma_PPY = c("ENSG00000108849"),
    ductal_KRT19 = c("ENSG00000171345"),
    acinar_PRSS1 = c("ENSG00000204983"),
    immune_CD45 = c("ENSG00000081237"),
    endothelial_VWF = c("ENSG00000110799"),
    unclear_SDS = c("ENSG00000135094"),
    stellate_RGS5 = c("ENSG00000232995", "ENSG00000143248")
  )
}

png(file.path(opts$outdir, "barplots_marker_genes_pre.png"))
plotDecontXMarkerPercentage(results,
    markers = markers,
    assayName = "counts")
dev.off()

png(file.path(opts$outdir, "barplots_marker_genes_post.png"))
plotDecontXMarkerPercentage(results,
    markers = markers,
    assayName = "decontXcounts")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Saving results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Saving complete decontX results.\n")
saveRDS(results, file=file.path(opts$outdir, "decontx_results.rds"))

cat("Saving updated count matrices.\n")
raw_counts_filtered = counts(x_nuclei)[,metadata(results)[[1]][["contamination"]] < opts$max_contamination]
decontaminated_counts_filtered = decontXcounts(results)[,metadata(results)[[1]][["contamination"]]< opts$max_contamination]
saveRDS(raw_counts_filtered, file=file.path(opts$outdir, "counts_low_contamination_raw.rds"))
saveRDS(decontaminated_counts_filtered, file=file.path(opts$outdir, "counts_low_contamination_decontaminated.rds"))

saveRDS(results, file=file.path(opts$outdir, "full_results.rds"))


write.table(contamination_estimates, file = file.path(opts$outdir, "contamination_estimates.tsv"), sep='\t', quote=FALSE, col.names=FALSE)

if (!is.na(opts$filter_log)) {
line = paste(opts$step_name, length(results$decontX_contamination), length(colnames(raw_counts_filtered)), sep=',')
write(line, opts$filter_log, append=TRUE)
}


cat("DONE.\n")
