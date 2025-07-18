# DecontX snakemake pipeline

Snakemake pipeline that corrects for ambient RNA expression in single-cell RNA sequencing (scRNA-seq) data using DecontX.

## Overview

This project contains the code to run a two-phase tissue-agnostic decontamination protocol, as well as configurations and wrapper scripts to run the pipeline on sample data. The pipeline supports multiple tissue types through configurable marker genes and parameters (mainly for visualization).

## Quick start
```{bash}
git clone https://github.com/grenkoca/snakemake-decontx.git
cd snakemake-decontx
conda env create --name sc_decontx -f environment.yaml # you may need to enter [Y]
conda activate sc_decontx
chmod +X ./run_pipeline.sh
./run_pipeline.sh --configfile config/pancreatic_islet_config.yaml # will run demo data
```

### Installation 
Either install the provided conda environment (`./environment.yaml`) or use the docker image (using `./run_docker.sh`. Requires pulling `letaylor/sc_decontx`).


### Input data
The expected input data is a folder containing standard 10x outputs. Each folder should contain the following standard folders:
- `[sample]/outs/filtered_feature_bc_matrix`
- `[sample]/outs/raw_feature_bc_matrix`
  
with each `[raw/filtered]_feature_bc_matrix` folder containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. For reference, see the [provided sample data](./data/).


### Configuration
To use your own data:
1) Copy the most relevant template configuration from `./config/*.yaml`. You will want to modify it to fit your project. The Following fields are most relevant: `name`, `input_dir_base`, `input_dir_format`, `tissue_type`, and `samples`. See below for descriptions.
2) Update `./config/marker_genes.yaml` to ensure that your tissues marker genes are correct. 
3) Run the pipeline using your new config file: `./run_pipeline.sh --configfile config/[name].yaml`


### Config fields:
- `name`: The run name. This will determine the output location.  
- `input_dir_base`: The parent directory of all your sample folders. (See provided example in `./data/`)
- `input_path_format`: The template path that the pipeline uses to find feature matrices. The first part should match `input_base_dir`, but with double brackets around the sample wildcard.
- `tissue_type`: A string corresponding to a tissue type in `./config/marker_genes.yaml`. Must match a key under `cell_type_markers`. This will define what cell types are present in your tissue + what genes are used as their markers. 
- `samples`: The names of samples. Should match folderes in the directory structure you have defined in `input_dir_base` and `input_path_format`


## Pipeline Workflow

The pipeline follows a step-by-step decontamination process with iterative clustering and decontamination rounds:

1. **prep_droplets**: Separates cell-containing droplets from empty droplets using the raw and filtered 10x data. Outputs nuclei counts and empty droplet counts.

2. **seurat_prelim**: Performs initial clustering on the nuclei data using Seurat with configurable resolution and PCA dimensions.

3. **decontx_prelim**: First round of decontamination using DecontX with the initial clusters. Removes highly contaminated cells based on the `max_contamination` threshold.

4. **seurat_round2**: Second round of clustering on the preliminary decontaminated data to refine cluster assignments.

5. **decontx_round2**: Second round of decontamination using the refined clusters. Generates final contamination estimates and produces the cleaned count matrix.

6. **seurat_round3**: Final clustering on the fully decontaminated data for downstream analysis.

7. **dump_seurat_object**: Exports the final cleaned data in standard 10x format (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).

## Configuration Parameters

### Decontamination Parameters
- `max_contamination`: Maximum contamination threshold (default: 0.1). Cells with contamination estimates above this value are removed.
- `delta_first`: First delta parameter for DecontX (default: 10). Controls the strength of the first decontamination round.
- `delta_second`: Second delta parameter for DecontX (default: 30). Controls the strength of the second decontamination round.

### Clustering Parameters
- `resolution`: Clustering resolution for Seurat (default: 0.8). Higher values result in more clusters.
- `pca_dims`: Number of PCA dimensions to use (default: 20).
- `min_pca_dim`: Minimum PCA dimension to include (default: 1).
- `max_pca_dim`: Maximum PCA dimension to include (default: 20).

### QC Metrics
The pipeline calculates quality control metrics for genes matching specific patterns:
- `mt`: Mitochondrial genes (pattern: `^MT-`)
- `ribo`: Ribosomal genes (pattern: `^RP[SL]`)
- Tissue-specific genes (e.g., for pancreatic islets: `INS` for insulin, `GCG` for glucagon)

### TODO:
Fix broken plots
- [ ] decontx: umap of pre/post marker expression are identical
- [ ] seurat: density plots for QC are empty
- [ ] seurat: combined feature umap only contains first 4 markers
- [ ] seurat: pca is colored by single value, `scRNA_decontamination`    

### Acknowledgements
This pipeline originally came from a collaboration between the labs of Dr. Francis Collins and Dr. Mike Boehnke. Special thanks to Drs. Cassie Robertson, Leland Taylor, Henry Taylor, and all others involved (please message me if you are left off)

