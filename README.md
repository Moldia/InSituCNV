# InSituCNV

Reusable notebook workflow for inferring copy-number variation (CNV) profiles
from image-based spatial transcriptomics data.

This repository contains the general package version of the workflow used in
***In Situ* Inference of Copy Number Variations in Image-Based Spatial
Transcriptomics** by Jensen et al. The manuscript reproduction code is kept in a
separate repository: https://github.com/Moldia/InSituCNV-manuscript

## Install

Clone the repository, create the conda environment, and install the package in
editable mode:

```bash
git clone https://github.com/Moldia/InSituCNV.git
cd InSituCNV
conda env create -f insitucnv.yml
conda activate insitucnv_env
pip install -e .
jupyter lab
```

Open the notebooks from the `notebooks/` directory in JupyterLab.

## What You Need Before Running

Prepare an `.h5ad` file with the information needed by the notebooks:

- raw counts, preferably in `adata.layers["raw_counts"]`; if raw counts are in
  `adata.X`, the first notebook can copy them into `adata.layers["raw_counts"]`;
- spatial coordinates in `adata.obsm["spatial"]`;
- a nearest-neighbor graph for smoothing, usually from `scanpy.pp.neighbors`;
- an `adata.obs` column that identifies normal or healthy reference cells for
  `infercnvpy`;
- the exact category names in that reference column that should be used as the
  normal reference, for example immune, stromal, or other non-tumor cell types;
- gene names that can be matched to genomic coordinates. The package includes
  the default `infercnvpy` gene coordinate table, and the notebooks show where
  to adjust this if your gene annotation differs.

Do not use normalized or log-transformed values as raw counts. The CNV workflow
normalizes, smooths, and log-transforms the raw counts itself.

## Notebook Workflow

Run the notebooks in order for your own dataset, editing the setup cells at the
top of each notebook.

### 1. Run InSituCNV

`notebooks/run_insitucnv.ipynb`

Use this notebook to run the full CNV analysis step by step:

- load your spatial transcriptomics `.h5ad` file;
- choose the raw count layer name with `RAW_LAYER`;
- choose the normal/reference annotation with `REFERENCE_KEY` and
  `REFERENCE_CATEGORIES`;
- normalize raw counts;
- smooth normalized counts over the neighbor graph;
- add genomic positions;
- run `infercnvpy`;
- cluster CNV profiles across selected Leiden resolutions;
- visualize chromosome heatmaps and spatial CNV cluster plots.

After reviewing the heatmap and spatial plot, manually edit:

- `NORMAL_CLUSTERS`;
- `TUMOR_CLUSTERS`, if you want to specify tumor clusters directly;
- `TUMOR_CLONE_CLUSTERS`, if some tumor clusters have distinct enough CNV
  profiles to report as separate tumor clones.

No cluster is marked as normal automatically.

## Outputs

The notebooks write results under their configured output directories, usually
`results/...` or `outputs/...`. Typical outputs include:

- checked input `.h5ad` files;
- AnnData objects with CNV values and CNV cluster labels;
- CNV chromosome heatmaps;
- spatial CNV cluster plots;
- optional manually annotated tumor/normal and tumor-clone cell group tables;
- optional mean CNV profiles for manually selected tumor cells.
