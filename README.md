# InSituCNV

Reusable notebook workflow for inferring copy-number variation (CNV) profiles
from image-based spatial transcriptomics data.

This repository contains the general package version of the workflow used in
***In Situ* Inference of Copy Number Variations in Image-Based Spatial
Transcriptomics** by Jensen et al. The manuscript reproduction code is kept in a
separate repository: https://github.com/Moldia/InSituCNV-manuscript

## References

InSituCNV builds upon several excellent packages in the single-cell and spatial
transcriptomics ecosystem:

- **[infercnvpy](https://github.com/icbi-lab/infercnvpy)**: Core CNV inference logic.
- **[scVelo](https://github.com/theislab/scvelo)**: Transcriptomic smoothing and analysis.
- **[Scanpy](https://github.com/scverse/scanpy)**: General single-cell analysis framework.
- **[AnnData](https://github.com/scverse/anndata)**: Data structures for single-cell data.

If you use InSituCNV in your research, please cite:

Jensen et al. *In Situ* inference of copy number variations in image-based spatial transcriptomics. bioRxiv (2025).
https://doi.org/10.1101/2025.07.02.662761

## Installation

```bash
pip install insitucnv        # once released on PyPI
```

From a clone (also see [CONTRIBUTING.md](CONTRIBUTING.md) for the full dev setup):

```bash
git clone https://github.com/Moldia/InSituCNV.git
cd InSituCNV
conda env create -f insitucnv.yml
conda activate insitucnv_env
```

## Quickstart

```python
import scanpy as sc
import insitucnv as icv

adata = sc.read_h5ad(icv.download_example_dataset())   # or your own .h5ad

result = icv.run_insitucnv(
    adata,
    output_dir="results/sample",
    reference_key="cell_type",
    reference_categories=["T_cells", "B_cells", "Stromal"],
)
```

Or work through `notebooks/run_insitucnv.ipynb`, which runs the same steps one at
a time and downloads the example dataset by default.

## What You Need Before Running

An `.h5ad` file with:

- **raw counts** in `adata.layers["raw_counts"]` (or in `adata.X` — the workflow
  copies them). Do not pass normalized or log-transformed values as raw counts;
  the workflow normalizes, smooths, and log-transforms them itself.
- **spatial coordinates** in `adata.obsm["spatial"]`;
- an **`adata.obs` column** identifying non-tumor reference cells (immune,
  stromal, …) and the category names to use, passed as `reference_key` and
  `reference_categories`;
- **human gene symbols** (GRCh38) in `adata.var_names`. Genomic coordinates come
  from the bundled `infercnvpy` table; pass `gene_reference_path=` a CSV with
  `gene_name,chromosome,start,end` for a different annotation.

A nearest-neighbor graph for smoothing is built automatically when absent
(`build_neighbors=True`); pass `build_neighbors=False` to require your own.

## Outputs and manual annotation

`run_insitucnv` writes the CNV `AnnData`, chromosome heatmaps, spatial cluster
plots and a `run_summary.json` under `output_dir`. Tumor/normal and clone labels
are **not** inferred automatically: after reviewing the heatmap, pass the cluster
labels you choose to `icv.tl.assign_cnv_status` /
`icv.tl.assign_cnv_subclones`, then export Xenium Explorer-compatible tables with
`icv.tl.export_cell_groups` and `icv.tl.export_mean_cnv_per_gene`.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, tests, linting,
docs and release process. Changes are recorded in [CHANGELOG.md](CHANGELOG.md).
