options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(ggplot2)
library(yaml)
theme_set(theme_bw(base_size = 12))

# Source configuration utilities
source("workflow/scripts/config_utils.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script reads in raw counts matrices for separately for each batch
# and generates preliminary cell type clusters using configurable parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

option_list <- list(
  make_option(
    c("--counts"), type = "character", help = "Input quality nuclei counts matrix."
  ),
  make_option(
    c("--resolution"), type = "numeric", help = "Clustering resolution."
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output destination for clustering results."
  ),
  make_option(
    c("--config"), type = "character", help = "Path to main configuration file."
  ),
  make_option(
    c("--gene_mapping"), type = "character", help = "Path to custom gene mapping RDS file."
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read configuration
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!is.null(opts$config) && file.exists(opts$config)) {
  main_config <- yaml.load_file(opts$config)
  
  # Read marker genes configuration
  marker_config_file <- main_config$marker_genes_config
  if (!file.exists(marker_config_file)) {
    marker_config_file <- file.path(dirname(opts$config), marker_config_file)
  }
  
  marker_config <- read_marker_genes_config(marker_config_file)
  qc_metrics <- get_qc_metrics(main_config)
  clustering_params <- get_clustering_params(main_config)
  tissue_type <- ifelse(is.null(main_config$tissue_type), "general", main_config$tissue_type)
} else {
  # Fallback to default parameters if no config provided
  qc_metrics <- list(
    list(name = "mt", pattern = "^MT-", pattern_type = "HGNC", description = "Mitochondrial genes"),
    list(name = "ribo", pattern = "^RP[SL]", pattern_type = "HGNC", description = "Ribosomal genes")
  )
  clustering_params <- list(
    resolution = 0.8,
    pca_dims = 20,
    min_pca_dim = 1,
    max_pca_dim = 20
  )
  tissue_type <- "general"
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
counts = readRDS(opts$counts)

### Create Seurat object
sobj <- CreateSeuratObject(counts = counts, project = "scRNA_decontamination", min.cells = 0, min.features = 0)
sobj

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC - Calculate configurable QC metrics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load custom gene mapping if provided
if (!is.null(opts$gene_mapping) && file.exists(opts$gene_mapping)) {
  cat("Loading custom gene mapping from:", opts$gene_mapping, "\n")
  gene_mapping <- readRDS(opts$gene_mapping)
  cat("Loaded", nrow(gene_mapping), "gene mappings\n")
} else {
  gene_mapping <- NULL
  cat("No custom gene mapping provided, using direct pattern matching\n")
}

for (metric in qc_metrics) {
  metric_name <- paste0("percent.", metric$name)
  
  # Check if pattern type is specified and is HGNC
  if (!is.null(metric$pattern_type) && metric$pattern_type == "HGNC") {
    # Use custom gene mapping for HGNC pattern matching
    if (!is.null(gene_mapping)) {
      ensembl_genes <- rownames(sobj)
      
      # Filter gene mapping to only include genes in our dataset
      available_mapping <- gene_mapping[gene_mapping$ensembl_id %in% ensembl_genes, ]
      
      # Find genes matching the pattern
      matching_symbols <- grep(metric$pattern, available_mapping$gene_symbol, value = TRUE)
      matching_ensembl <- available_mapping$ensembl_id[available_mapping$gene_symbol %in% matching_symbols]
      
      # Calculate percentage using ENSGID names
      if (length(matching_ensembl) > 0) {
        sobj[[metric_name]] <- PercentageFeatureSet(sobj, features = matching_ensembl)
        cat("Found", length(matching_ensembl), "genes matching pattern", metric$pattern, "for metric", metric$name, "\n")
      } else {
        warning(paste("No genes found matching pattern", metric$pattern, "for metric", metric$name))
        sobj[[metric_name]] <- 0
      }
    } 
  } else {
    # Use pattern directly (for ENSGID patterns or legacy configs)
    sobj[[metric_name]] <- PercentageFeatureSet(sobj, pattern = metric$pattern)
  }
}

# Output QC plots for each configured metric
for (metric in qc_metrics) {
  metric_name <- paste0("percent.", metric$name)
  plot_file <- file.path(opts$outdir, paste0("seurat_density_", metric$name, ".png"))
  
  png(plot_file)
  p <- ggplot(sobj@meta.data) + 
    geom_density(aes_string(x = paste0(metric_name, "+0.01"))) + 
    scale_x_continuous(trans='log10') + 
    xlab(paste("Percent", metric$description)) +
    ggtitle(paste("Density plot for", metric$description))
  print(p)
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normalization and scaling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sobj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file.path(opts$outdir,"seurat_variable_features.png"), width=800, height=400)
plot1 + plot2
dev.off()

all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))

# Print the first 5 principal components
print(sobj[["pca"]], dims = 1:5, nfeatures = 5)

png(file.path(opts$outdir,"seurat_pca.png"), width=800, height=600)
VizDimLoadings(sobj, dims = 1:2, reduction = "pca")
dev.off()

png(file.path(opts$outdir,"seurat_pca_plot.png"), width=600, height=400)
DimPlot(sobj, reduction = "pca")
dev.off()

png(file.path(opts$outdir,"seurat_pca_heatmap.png"), width=800, height=600)
DimHeatmap(sobj, dims = 1, cells = 500, balanced = TRUE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Determine dimensionality
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
png(file.path(opts$outdir,"seurat_elbow_plot.png"), width=600, height=400)
ElbowPlot(sobj)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clustering nuclei using configurable parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dim_range <- clustering_params$min_pca_dim:min(clustering_params$max_pca_dim, clustering_params$pca_dims)
sobj <- FindNeighbors(sobj, dims = dim_range)
sobj <- FindClusters(sobj, resolution = opts$resolution)

# View clusters for first 5 nuclei
head(Idents(sobj), 5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# UMAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sobj <- RunUMAP(sobj, dims = dim_range)

png(file.path(opts$outdir,"seurat_umap.png"), width=600, height=400)
DimPlot(sobj, reduction = "umap")
dev.off()

# Create UMAP plot colored by QC metrics
for (metric in qc_metrics) {
  metric_name <- paste0("percent.", metric$name)
  plot_file <- file.path(opts$outdir, paste0("seurat_umap_", metric$name, ".png"))
  
  png(plot_file, width=600, height=400)
  p <- FeaturePlot(sobj, features = metric_name) + 
    ggtitle(paste("UMAP colored by", metric$description))
  print(p)
  dev.off()
}

# Save Seurat object
saveRDS(sobj, file = file.path(opts$outdir,"seurat_obj.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save clusters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d = sobj[[]]
d$barcode_gex = row.names(d)
d$clusters = paste0("Seurat_",d$seurat_clusters)
write.csv(d[,c("barcode_gex", "clusters")], file=file.path(opts$outdir, "seurat_clusters.csv"), row.names=FALSE, quote=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Feature plots for known markers (if available)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (exists("marker_config") && !is.null(marker_config)) {
  known_markers <- get_known_markers(marker_config, tissue_type)
  
  # Check which markers are present in the data
  available_markers <- intersect(known_markers, rownames(sobj))
  
  if (length(available_markers) > 0) {
    # Create feature plots for available markers
    for (marker in available_markers) {
      plot_file <- file.path(opts$outdir, paste0("seurat_feature_", marker, ".png"))
      
      png(plot_file, width=600, height=400)
      p <- FeaturePlot(sobj, features = marker) + 
        ggtitle(paste("Expression of", marker))
      print(p)
      dev.off()
    }
    
    # Create a combined feature plot for top markers

    width <- ceiling(sqrt(length(available_markers))) * 300 + 200
    height <- floor(sqrt(length(available_markers))) * 300
    png(file.path(opts$outdir, "seurat_feature_combined.png"), width=width, height=height)
    p <- FeaturePlot(sobj, features = available_markers)
    print(p)
    dev.off()
  }
}

cat("Seurat analysis completed successfully!\n")
cat("Output files saved to:", opts$outdir, "\n")
cat("Clustering resolution used:", opts$resolution, "\n")
cat("PCA dimensions used:", paste(dim_range, collapse=":"), "\n")
