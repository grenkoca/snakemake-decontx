#!/usr/bin/env Rscript

# Script to build comprehensive ENSGID to HGNC mapping from 10x features.tsv.gz files
# This creates a reliable mapping table directly from the input data

library(optparse)
library(dplyr)
library(data.table)

option_list <- list(
  make_option(
    c("--input_dirs"), type = "character", 
    help = "Comma-separated list of input directories containing features.tsv.gz files"
  ),
  make_option(
    c("--output"), type = "character", 
    help = "Output file path for the gene mapping RDS file"
  ),
  make_option(
    c("--report"), type = "character", 
    help = "Output file path for the mapping report"
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# Parse input directories
input_dirs <- strsplit(opts$input_dirs, ",")[[1]]
input_dirs <- trimws(input_dirs)

cat("Building gene mapping from", length(input_dirs), "input directories\n")

# Initialize data structures
all_mappings <- data.table()
file_info <- data.table()

# Process each input directory
for (i in seq_along(input_dirs)) {
  dir_path <- input_dirs[i]
  features_file <- file.path(dir_path, "features.tsv.gz")
  
  if (!file.exists(features_file)) {
    cat("Warning: features.tsv.gz not found in", dir_path, "\n")
    next
  }
  
  cat("Processing", features_file, "\n")
  
  # Read features file
  features <- fread(features_file, header = FALSE, 
                    col.names = c("ensembl_id", "gene_symbol", "feature_type"))
  
  # Filter to Gene Expression features only
  features <- features[feature_type == "Gene Expression"]
  
  # Add source information
  features$source_dir <- dir_path
  features$source_file <- i
  
  # Store file info
  file_info <- rbind(file_info, data.table(
    file_index = i,
    directory = dir_path,
    n_genes = nrow(features),
    n_unique_ensembl = length(unique(features$ensembl_id)),
    n_unique_symbols = length(unique(features$gene_symbol))
  ))
  
  # Combine with all mappings
  all_mappings <- rbind(all_mappings, features)
}

cat("Total mappings collected:", nrow(all_mappings), "\n")

# Create comprehensive mapping table
gene_mapping <- all_mappings[, .(
  gene_symbol = first(gene_symbol),
  n_files = .N,
  n_unique_symbols = length(unique(gene_symbol)),
  all_symbols = paste(unique(gene_symbol), collapse = "|"),
  source_files = paste(unique(source_file), collapse = ",")
), by = ensembl_id]

# Add conflict flags
gene_mapping[, symbol_conflict := n_unique_symbols > 1]
gene_mapping[, multiple_files := n_files > 1]

# Check for reverse conflicts (same symbol mapping to multiple ENSGs)
symbol_conflicts <- gene_mapping[, .(
  n_ensembl = .N,
  n_unique_ensembl = length(unique(ensembl_id)),
  all_ensembl = paste(unique(ensembl_id), collapse = "|")
), by = gene_symbol]

symbol_conflicts[, ensembl_conflict := n_unique_ensembl > 1]

# Merge back symbol conflict info
gene_mapping <- merge(gene_mapping, 
                     symbol_conflicts[, .(gene_symbol, ensembl_conflict)], 
                     by = "gene_symbol", all.x = TRUE)

# Summary statistics
n_total_genes <- nrow(gene_mapping)
n_symbol_conflicts <- sum(gene_mapping$symbol_conflict)
n_ensembl_conflicts <- sum(gene_mapping$ensembl_conflict, na.rm = TRUE)
n_clean_mappings <- sum(!gene_mapping$symbol_conflict & !gene_mapping$ensembl_conflict)

cat("\n=== GENE MAPPING SUMMARY ===\n")
cat("Total unique ENSEMBL IDs:", n_total_genes, "\n")
cat("Clean 1:1 mappings:", n_clean_mappings, "\n")
cat("ENSEMBL IDs with multiple symbols:", n_symbol_conflicts, "\n")
cat("Gene symbols with multiple ENSEMBL IDs:", n_ensembl_conflicts, "\n")
cat("Mapping quality:", round(n_clean_mappings / n_total_genes * 100, 1), "%\n")

# Create final mapping table for use in analysis
# For conflicts, we'll use the most common symbol/ENSEMBL pairing
final_mapping <- gene_mapping[, .(
  ensembl_id = ensembl_id,
  gene_symbol = gene_symbol,
  has_conflicts = symbol_conflict | ensembl_conflict,
  n_files = n_files
)]

# Save the mapping
saveRDS(final_mapping, file = opts$output)
cat("Gene mapping saved to:", opts$output, "\n")

# Generate detailed report
if (!is.null(opts$report)) {
  
  # Create detailed report
  report_lines <- c(
    "=== GENE MAPPING REPORT ===",
    paste("Generated on:", Sys.time()),
    paste("Input directories:", length(input_dirs)),
    "",
    "=== FILE SUMMARY ===",
    paste(sprintf("%-30s %8s %12s %12s", "Directory", "N_genes", "N_ensembl", "N_symbols")),
    paste(sprintf("%-30s %8d %12d %12d", 
                  basename(file_info$directory), 
                  file_info$n_genes, 
                  file_info$n_unique_ensembl, 
                  file_info$n_unique_symbols)),
    "",
    "=== MAPPING SUMMARY ===",
    paste("Total unique ENSEMBL IDs:", n_total_genes),
    paste("Clean 1:1 mappings:", n_clean_mappings),
    paste("ENSEMBL IDs with multiple symbols:", n_symbol_conflicts),
    paste("Gene symbols with multiple ENSEMBL IDs:", n_ensembl_conflicts),
    paste("Mapping quality:", paste0(round(n_clean_mappings / n_total_genes * 100, 1), "%")),
    ""
  )
  
  # Add examples of conflicts if any
  if (n_symbol_conflicts > 0) {
    report_lines <- c(report_lines,
      "=== SYMBOL CONFLICTS (ENSEMBL -> multiple symbols) ===",
      paste("Top 10 examples:"),
      paste(sprintf("%-20s -> %s", 
                    gene_mapping[symbol_conflict == TRUE][1:min(10, .N)]$ensembl_id,
                    gene_mapping[symbol_conflict == TRUE][1:min(10, .N)]$all_symbols)),
      ""
    )
  }
  
  if (n_ensembl_conflicts > 0) {
    report_lines <- c(report_lines,
      "=== ENSEMBL CONFLICTS (symbol -> multiple ENSEMBL) ===",
      paste("Top 10 examples:"),
      paste(sprintf("%-20s -> %s", 
                    symbol_conflicts[ensembl_conflict == TRUE][1:min(10, .N)]$gene_symbol,
                    symbol_conflicts[ensembl_conflict == TRUE][1:min(10, .N)]$all_ensembl)),
      ""
    )
  }
  
  # Check for MT and ribosomal genes
  mt_genes <- final_mapping[grepl("^MT-", gene_symbol)]
  ribo_genes <- final_mapping[grepl("^RP[SL]", gene_symbol)]
  
  report_lines <- c(report_lines,
    "=== QC GENE DETECTION ===",
    paste("Mitochondrial genes (MT-*):", nrow(mt_genes)),
    paste("Ribosomal genes (RP[SL]*):", nrow(ribo_genes)),
    ""
  )
  
  if (nrow(mt_genes) > 0) {
    report_lines <- c(report_lines,
      "Sample mitochondrial genes:",
      paste(sprintf("  %s -> %s", mt_genes$ensembl_id[1:min(5, nrow(mt_genes))], 
                                   mt_genes$gene_symbol[1:min(5, nrow(mt_genes))]))
    )
  }
  
  if (nrow(ribo_genes) > 0) {
    report_lines <- c(report_lines,
      "Sample ribosomal genes:",
      paste(sprintf("  %s -> %s", ribo_genes$ensembl_id[1:min(5, nrow(ribo_genes))], 
                                   ribo_genes$gene_symbol[1:min(5, nrow(ribo_genes))]))
    )
  }
  
  # Write report
  writeLines(report_lines, opts$report)
  cat("Detailed report saved to:", opts$report, "\n")
}

cat("Gene mapping build completed successfully!\n")