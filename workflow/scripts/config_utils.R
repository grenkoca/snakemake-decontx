# Configuration utilities for tissue-agnostic scRNA-seq decontamination pipeline
# This script provides functions to read and parse configuration files

library(yaml)

# Function to read marker genes configuration
read_marker_genes_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop(paste("Marker genes configuration file not found:", config_file))
  }
  
  config <- yaml.load_file(config_file)
  return(config)
}

# Function to get known markers for a specific tissue type
get_known_markers <- function(marker_config, tissue_type = "general") {
  
  # Get general known markers
  known_markers <- marker_config$known_markers
  
  # Add tissue-specific markers if available
  if (tissue_type %in% names(marker_config$cell_type_markers)) {
    tissue_markers <- marker_config$cell_type_markers[[tissue_type]]
    
    # Extract all markers from all cell types for this tissue
    tissue_marker_list <- unlist(tissue_markers, use.names = FALSE)
    known_markers <- c(known_markers, tissue_marker_list)
  }
  
  # Remove duplicates
  known_markers <- unique(known_markers)
  
  return(known_markers)
}

# Function to get cell type markers for a specific tissue type
get_cell_type_markers <- function(marker_config, tissue_type = "general") {
  
  if (tissue_type %in% names(marker_config$cell_type_markers)) {
    return(marker_config$cell_type_markers[[tissue_type]])
  } else if (tissue_type == "general") {
    return(marker_config$general)
  } else {
    warning(paste("Tissue type", tissue_type, "not found in marker configuration. Using general markers."))
    return(marker_config$general)
  }
}

# Function to create markers list in the format expected by the decontamination scripts
create_markers_list <- function(cell_type_markers) {
  markers_list <- list()
  
  for (cell_type in names(cell_type_markers)) {
    # Create a name that follows the pattern: celltype_marker
    marker_name <- paste0(cell_type, "_markers")
    markers_list[[marker_name]] <- cell_type_markers[[cell_type]]
  }
  
  return(markers_list)
}

# Function to read QC metrics configuration
get_qc_metrics <- function(config) {
  if ("qc_metrics" %in% names(config)) {
    return(config$qc_metrics)
  } else {
    # Default QC metrics if not specified
    default_qc <- list(
      list(name = "mt", pattern = "^MT-", description = "Mitochondrial genes"),
      list(name = "ribo", pattern = "^RP[SL]", description = "Ribosomal genes")
    )
    return(default_qc)
  }
}

# Function to get clustering parameters
get_clustering_params <- function(config) {
  if ("clustering" %in% names(config)) {
    return(config$clustering)
  } else {
    # Default clustering parameters
    default_clustering <- list(
      resolution = 0.8,
      pca_dims = 20,
      min_pca_dim = 1,
      max_pca_dim = 20
    )
    return(default_clustering)
  }
}