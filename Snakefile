#!/usr/bin/env python3

from os.path import join
import os
import pandas as pd
from functools import partial
import math
import glob
from pathlib import Path

name = config["name"]
_data = partial(os.path.join, "data")
_results = partial(os.path.join, "results", name)
_resources = partial(os.path.join, "resources")
_logs = partial(_results, "logs")

### Samples to process (fiveprime.yaml)
samples = config["samples"]
wildcard_constraints:
   sample="|".join(samples)

# Added by Caleb (02-02-2023)
INPUT_PATH_FORMAT = config["input_path_format"]
INPUT_DIR_BASE    = config["input_dir_base"]
INPUT_RAW      = Path(INPUT_DIR_BASE) / "raw_feature_bc_matrix"
INPUT_FILTERED = Path(INPUT_DIR_BASE) / "filtered_feature_bc_matrix"


INPUT_MTX      = Path(INPUT_PATH_FORMAT.format(data_file = 'matrix.mtx.gz'))
INPUT_BARCODES = Path(INPUT_PATH_FORMAT.format(data_file = 'barcodes.tsv.gz'))
INPUT_FEATURES = Path(INPUT_PATH_FORMAT.format(data_file = 'features.tsv.gz'))


rule all:
    input:
        expand(_results("cleaned/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz"), sample=samples),
        expand(_results("cleaned/{sample}/outs/raw_feature_bc_matrix/features.tsv.gz"), sample=samples),
        expand(_results("cleaned/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz"), sample=samples)


# Figure out which droplets contain cells and which droplets contain nuclei
rule prep_droplets:
    input:
        raw_data = directory(INPUT_RAW),
        filtered_data = directory(INPUT_FILTERED)
    output:
        counts_nuclei = _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_protein_coding/{sample}/counts_empty.rds"),
        filter_file = _results("filter_counts/{sample}/cells_filtered.csv")
    params:
        outdir = lambda wildcards: _results(f"counts_protein_coding/{wildcards.sample}/"),
        outdir2 = lambda wildcards: _results(f"filter_counts/{wildcards.sample}"),
    shell:
        "mkdir -p {params.outdir}; "
        "mkdir -p {params.outdir2}; "
        "echo step,start,end > {params.outdir2}/cells_filtered.csv; "
        """
        {config[Rscript_binary]} workflow/scripts/get_10x_empty_droplets.R \
            --raw_10x_dir {input.raw_data} \
            --filtered_10x_dir {input.filtered_data} \
            --filter_log {output.filter_file} \
            --counts_nuclei {output.counts_nuclei} \
            --counts_empty {output.counts_empty}
        """

rule seurat_prelim:
    input:
        counts = _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
        original_features = INPUT_FEATURES,
        original_barcodes = INPUT_BARCODES
    output:
        #_results("seurat_prelim/{sample}/seurat_obj.rds"),
        _results("seurat_prelim/{sample}/seurat_clusters.csv"),
    params:
        outdir = lambda wildcards:  _results(f"seurat_prelim/{wildcards.sample}"),
        resolution = 0.8
    shell:
        "mkdir -p {params.outdir}; "
        """
        {config[Rscript_binary]} workflow/scripts/run_seurat.R.bak \
            --counts {input.counts} \
            --resolution {params.resolution} \
            --outdir {params.outdir} \
        """


rule decontx_prelim:
    input:
        counts_nuclei = _results("counts_protein_coding/{sample}/counts_nuclei.rds"),
        counts_empty = _results("counts_protein_coding/{sample}/counts_empty.rds"),
        clusters = _results("seurat_prelim/{sample}/seurat_clusters.csv"),
    output:
        _results("decontx_prelim/{sample}/counts_low_contamination_raw.rds"),
    params:
        sample = lambda wildcards: wildcards.sample,
        outdir = lambda wildcards: _results(f"decontx_prelim/{wildcards.sample}"),
        max_contamination = config["max_contamination"],
        delta_first = 10,
        delta_second = 30
    shell:
        "mkdir -p {params.outdir}; "
        """
        {config[Rscript_binary]} workflow/scripts/run_decontx.R \
            --counts_nuclei {input.counts_nuclei} \
            --counts_empty {input.counts_empty} \
            --clusters {input.clusters} \
            --max_contamination {params.max_contamination} \
            --delta_first {params.delta_first} \
            --delta_second {params.delta_second} \
            --outdir {params.outdir}
        """

rule seurat_round2:
    input:
        counts = _results("decontx_prelim/{sample}/counts_low_contamination_raw.rds"),
        original_features = INPUT_FEATURES,
    output:
        #_results("seurat_round2/{sample}/seurat_obj.rds"),
        _results("seurat_round2/{sample}/seurat_clusters.csv"),
    params:
        outdir = lambda wildcards: _results(f"seurat_round2/{wildcards.sample}"),
        resolution = 0.8,
    shell:
        "mkdir -p {params.outdir}; "
        """
        {config[Rscript_binary]} workflow/scripts/run_seurat.R.bak \
            --counts {input.counts} \
            --resolution {params.resolution} \
            --outdir {params.outdir} \
        """


rule decontx_round2:
    input:
        counts_nuclei = _results("decontx_prelim/{sample}/counts_low_contamination_raw.rds"),
        counts_empty = _results("counts_protein_coding/{sample}/counts_empty.rds"),
        clusters = _results("seurat_round2/{sample}/seurat_clusters.csv"),
    output:
        results = _results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"),
        contamination  = _results("decontx_round2/{sample}/contamination_estimates.tsv")
    conda:
        "Renv"
    params:
        sample = lambda wildcards: wildcards.sample,
        outdir = lambda wildcards: _results(f"decontx_round2/{wildcards.sample}"),
        max_contamination = config["max_contamination"],
        delta_first = 10,
        delta_second = 30
    shell:
        "mkdir -p {params.outdir}; "
        """
        {config[Rscript_binary]} workflow/scripts/run_decontx.R \
            --counts_nuclei {input.counts_nuclei} \
            --counts_empty {input.counts_empty} \
            --clusters {input.clusters} \
            --max_contamination {params.max_contamination} \
            --delta_first {params.delta_first} \
            --delta_second {params.delta_second} \
            --outdir {params.outdir}; 
        """
        "{config[python_binary]} workflow/scripts/modify_df.py"
            " --input {output.contamination}"
            " --output {output.contamination}"
            " --new_columns experiment_id"
            " --new_values {params.sample}; " 

rule seurat_round3:
    input:
        counts = _results("decontx_round2/{sample}/counts_low_contamination_decontaminated.rds"),
        original_features = INPUT_FEATURES,
    output:
        _results("seurat_round3/{sample}/seurat_obj.rds"),
        _results("seurat_round3/{sample}/seurat_clusters.csv"),
        directory(_results("seurat_round3/{sample}"))
    params:
        outdir = lambda wildcards: _results(f"seurat_round3/{wildcards.sample}"),
        resolution = 0.8,
    shell:
        "mkdir -p {params.outdir}; "
        """
        {config[Rscript_binary]} workflow/scripts/run_seurat.R.bak \
            --counts {input.counts} \
            --resolution {params.resolution} \
            --outdir {params.outdir} \
        """


checkpoint dump_seurat_object:
    input:
        rds = _results("seurat_round3/{sample}/seurat_obj.rds"),
        old = INPUT_FEATURES
    output:
        _results("cleaned/{sample}/outs/raw_feature_bc_matrix/barcodes.tsv.gz"),
        _results("cleaned/{sample}/outs/raw_feature_bc_matrix/features.tsv.gz"),
        _results("cleaned/{sample}/outs/raw_feature_bc_matrix/matrix.mtx.gz"),
        directory(_results("cleaned/{sample}")),
        directory(_results("cleaned/{sample}/outs/raw_feature_bc_matrix/"))
    params:
        outdir = lambda wildcards: _results(f"cleaned/{wildcards.sample}/outs/raw_feature_bc_matrix"),
        outdir_parent = lambda wildcards: _results(f"cleaned/{wildcards.sample}/outs/")
    shell:
        "rm -rfv {params.outdir}; "
        "mkdir -p {params.outdir_parent}; "
        "{config[Rscript_binary]} workflow/scripts/dump_seurat_RDS.R"
            " {input.rds} "
            " {params.outdir}; "
        #"gzip {params.outdir}/*; "
        "{config[Rscript_binary]} workflow/scripts/update_ensgid.R"
            " {params.outdir}/features.tsv.gz"
            " {input.old}"
            " {params.outdir}/features.tsv.gz; "


