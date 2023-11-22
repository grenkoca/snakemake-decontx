# DecontX snakemake pipeline

Snakemake pipeline that corrects for ambient RNA expression in single-cell RNA sequencing (scRNA-seq) data using DecontX.

## Overview

This project contains the code to run a two-phase islet decontamination protocol, as well as configurations and wrapper scripts to run the pipeline on sample data.

### Quickstart using Docker

1. Clone the GitHub repo and cd into the repo directory.

```bash
# clone the repo
git clone https://github.com/CollinsLabBioComp/islet_decontamination.git

# set code base path
SNK_REPO="$(pwd)/islet_decontamination"
cd ${SNK_REPO}
```

2. Launch the Docker app and download the Docker image
```bash
docker pull letaylor/sc_decontx:latest
```

3. Run the sample data:
```
chmod +x ./run_docker.sh
./run_docker.sh
```

### Input data
The expected input data is a folder containing standard 10x outputs. Each folder should contain the following standard folders:
- `[sample]/outs/filtered_feature_bc_matrix`
- `[sample]/outs/raw_feature_bc_matrix`
  
with each `[raw/filtered]_feature_bc_matrix` folder containing `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`. For reference, see the [provided sample data](https://github.com/CollinsLabBioComp/islet_decontamination/tree/main/data).


### Configuration
To configure to use your own data:

1) Update `workflow/src/threeprime.yaml` to change the run ID (`name`) and specify the sample IDs (`samples`). _note: if your 10x output directory format differs, you may need to update `input_dir_basename` and `input_path_format` to match._
2) Place the 10x output folder (containing a minimum of `outs/`, see [here](https://github.com/CollinsLabBioComp/islet_decontamination/tree/main/data) for reference) in the `./data/` folder. Alternatively, you can modify `input_dir_base` and `input_path_format` so the base (_i.e._ `data/`) points to your parent directory containing all samples.
