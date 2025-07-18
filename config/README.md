# Configuration Guide for scRNA-seq Decontamination Pipeline

This directory contains configuration files for running the decontamination pipeline with different tissue types.

## Available Configurations

### Pre-configured Tissue Types

- **`pancreatic_islet_config.yaml`** - Pancreatic islet-specific configuration (original)
- **`brain_config.yaml`** - Brain/neural tissue configuration
- **`liver_config.yaml`** - Liver tissue configuration
- **`lung_config.yaml`** - Lung tissue configuration
- **`general_config.yaml`** - General tissue-agnostic configuration

### Usage

To use a specific configuration, run snakemake with the `--configfile` parameter:

```bash
snakemake --configfile config/brain_config.yaml
```

## Configuration Parameters

### Core Parameters

- **`name`** - Project name for outputs
- **`tissue_type`** - Tissue type identifier (used for marker gene selection)
- **`samples`** - List of sample IDs to process

### Decontamination Parameters

- **`max_contamination`** - Maximum contamination fraction threshold (0.1 = 10%)
- **`delta_first`** - First delta parameter for decontX
- **`delta_second`** - Second delta parameter for decontX

### Clustering Parameters

- **`resolution`** - Clustering resolution (higher = more clusters)
- **`pca_dims`** - Number of principal components to compute
- **`min_pca_dim`** - Minimum PC dimension to use
- **`max_pca_dim`** - Maximum PC dimension to use

### QC Metrics

Define custom QC metrics to calculate:

```yaml
qc_metrics:
  - name: "mt"
    pattern: "^MT-"
    description: "Mitochondrial genes"
  - name: "custom_gene"
    pattern: "^GENE_NAME$"
    description: "Custom gene description"
```

## Marker Genes Configuration

The `marker_genes.yaml` file defines tissue-specific marker genes for visualization and analysis.

### Structure

```yaml
# General markers used across all tissues
known_markers:
  - "ENSG00000081237"  # PTPRC (CD45)

# Tissue-specific markers
cell_type_markers:
  tissue_name:
    cell_type1:
      - "ENSG00000123456"  # Gene1
    cell_type2:
      - "ENSG00000789012"  # Gene2
```

### Adding New Tissue Types

1. Add your tissue type to the `cell_type_markers` section in `marker_genes.yaml`
2. Define cell type-specific marker genes using Ensembl IDs
3. Create a new configuration file (e.g., `my_tissue_config.yaml`)
4. Set `tissue_type: "my_tissue"` in your configuration

## Example: Creating a Kidney Configuration

1. **Add markers to `marker_genes.yaml`:**

```yaml
cell_type_markers:
  kidney:
    podocytes:
      - "ENSG00000115414"  # NPHS1 (Nephrin)
    tubular_cells:
      - "ENSG00000169174"  # PCSK1
    mesangial_cells:
      - "ENSG00000188404"  # SELL
```

2. **Create `kidney_config.yaml`:**

```yaml
name: "kidney_decontamination"
tissue_type: "kidney"
max_contamination: 0.1
clustering:
  resolution: 0.7
  pca_dims: 25
qc_metrics:
  - name: "mt"
    pattern: "^MT-"
    description: "Mitochondrial genes"
  - name: "NPHS1"
    pattern: "^NPHS1$"
    description: "Nephrin podocyte marker"
marker_genes_config: "config/marker_genes.yaml"
samples:
  - sample1
  - sample2
```

3. **Run the pipeline:**

```bash
snakemake --configfile config/kidney_config.yaml
```

## Tips for Parameter Optimization

### Clustering Resolution
- **Lower resolution (0.3-0.6)** - Fewer, broader clusters (good for well-defined cell types)
- **Higher resolution (0.8-1.5)** - More, finer clusters (good for discovering rare cell types)

### Contamination Thresholds
- **Lower threshold (0.05-0.1)** - More stringent filtering
- **Higher threshold (0.15-0.2)** - More permissive filtering

### Delta Parameters
- **Lower delta values** - More conservative decontamination
- **Higher delta values** - More aggressive decontamination

## Troubleshooting

1. **No markers found**: Ensure Ensembl IDs are correct and present in your data
2. **Poor clustering**: Adjust resolution or PCA dimensions
3. **Over-filtering**: Increase `max_contamination` threshold
4. **Under-filtering**: Decrease `max_contamination` threshold or adjust delta parameters