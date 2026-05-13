# InSituCNV

Reusable Python tools for inferring copy-number variation (CNV) profiles from
image-based spatial transcriptomics data.

This repository contains the general package version of the workflow used in
***In Situ* Inference of Copy Number Variations in Image-Based Spatial
Transcriptomics** by Jensen et al. The manuscript reproduction code is kept in a
separate repository: https://github.com/Moldia/InSituCNV-manuscript

## What The Package Does

`insitucnv` wraps the repeated notebook workflow into reusable functions:

- normalize raw counts;
- smooth normalized counts over a nearest-neighbor graph;
- normalize and log-transform the smoothed matrix;
- add genomic coordinates required by `infercnvpy`;
- run `infercnvpy`;
- cluster CNV profiles across one or more Leiden resolutions;
- plot CNV heatmaps and spatial CNV cluster maps;
- annotate CNV clusters as tumor/normal or subclones;
- export mean CNV profiles per gene.

## Installation

```bash
git clone https://github.com/Moldia/InSituCNV.git
cd InSituCNV
conda env create -f insitucnv.yml
conda activate insitucnv_env
pip install -e .
```

For an existing environment:

```bash
pip install -e .
```

## Minimal Python Workflow

For a full worked example on the 1105_BL Xenium breast cancer sample, see
`notebooks/03_tutorial_1105_BL_InSituCNV.ipynb`.

```python
import scanpy as sc
import insitucnv as icv

adata = sc.read_h5ad("sample_preprocessed.h5ad")

# Required inputs:
# - raw counts in adata.layers["raw_counts"] or adata.X
# - reference cell labels in adata.obs["cell_type"]
# - a neighbor graph for smoothing, for example from sc.pp.neighbors(...)
# - spatial coordinates in adata.obsm["spatial"] for spatial plots

adata = icv.tl.prepare_cnv_input(
    adata,
    raw_layer="raw_counts",
    target_sum=1e4,
    smoothing_neighbors=100,
)

icv.tl.run_infercnv(
    adata,
    reference_key="cell_type",
    reference_categories=["T_cells", "B_cells", "Myeloid", "Plasma", "Fibroblast", "Endothelial"],
    window_size=60,
    step=10,
    lfc_clip=4,
)

icv.tl.compute_cnv_neighbors(adata)
cluster_keys = icv.tl.cluster_cnv_resolutions(adata, resolutions=[0.1, 0.2, 0.3])

primary_key = "cnv_leiden_res0.1"
icv.tl.assign_cnv_status(adata, primary_key)

icv.pl.plot_chromosome_heatmap(adata, groupby=primary_key, output_path="plots/cnv_heatmap.png")
icv.pl.plot_spatial(adata, color=primary_key, output_path="plots/spatial_cnv_clusters.png")
icv.pl.plot_spatial(adata, color="cnv_status", output_path="plots/spatial_cnv_status.png")

icv.tl.export_mean_cnv_per_gene(adata, "results/tumor_mean_cnv_per_gene.tsv")
adata.write("results/sample_cnv.h5ad", compression="gzip")
```

## One-Command Workflow

For a prepared `.h5ad`:

```bash
insitucnv \
  --input-h5ad sample_preprocessed.h5ad \
  --sample-id sample_01 \
  --output-dir results/sample_01 \
  --reference-key cell_type \
  --reference-categories T_cells,B_cells,Myeloid,Plasma,Fibroblast,Endothelial \
  --smooth-neighbors 100 \
  --cluster-resolutions 0.1,0.2,0.3
```

For a Xenium output directory:

```bash
insitucnv \
  --input-dir /path/to/xenium/output \
  --sample-id sample_01 \
  --output-dir results/sample_01 \
  --annotation-csv cell_type_annotations.csv \
  --reference-key cell_type
```

## Outputs

The high-level workflow writes:

- `adata_cnv.h5ad`: AnnData with CNV results and cluster annotations;
- `run_summary.json`: parameters and selected keys;
- `cluster_resolution_metrics.csv`: optional metrics when requested;
- `tumor_mean_cnv_per_gene.tsv`: mean CNV signal for tumor cells;
- `plots/*_heatmap.png`: chromosome heatmaps;
- `plots/*_spatial.png`: spatial maps for each CNV clustering resolution.

## Notes

`infercnvpy` requires gene coordinates in `adata.var["chromosome"]`,
`adata.var["start"]`, and `adata.var["end"]`. By default InSituCNV uses the
same `infercnvpy` reference gene table used in the manuscript notebooks. You can
also pass your own gene annotation table through
`icv.pp.add_genomic_positions(..., reference_path="genes.tsv")`.
