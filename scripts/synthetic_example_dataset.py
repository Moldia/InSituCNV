"""Write a small synthetic spatial dataset for CI / offline notebook runs.

    python scripts/synthetic_example_dataset.py OUT.h5ad

The output has raw counts, ``obsm['spatial']`` and a ``cell_type`` column with
reference (T_cells/B_cells/Stromal) and Epithelial cells that carry a chr7 gain
and a chr10 loss, so the CNV workflow produces a visible signal. This is *not*
the published example dataset (see make_example_dataset.py) - it is a stand-in so
the notebook can run without network access.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def build() -> ad.AnnData:
    rng = np.random.default_rng(0)
    n_genes_per_chrom = 60
    chroms = [f"chr{i}" for i in range(1, 11)]
    genes, gene_chrom = [], []
    for chrom in chroms:
        for i in range(n_genes_per_chrom):
            genes.append(f"{chrom}_g{i}")
            gene_chrom.append(chrom)
    gene_chrom = np.array(gene_chrom)

    n_ref = {"T_cells": 300, "B_cells": 180, "Stromal": 140}
    n_epi = 500
    labels = sum(([k] * v for k, v in n_ref.items()), []) + ["Epithelial"] * n_epi
    n_obs = len(labels)

    counts = rng.poisson(1.5, size=(n_obs, len(genes))).astype(np.float32)
    epi = np.array(labels) == "Epithelial"
    counts[np.ix_(epi, gene_chrom == "chr7")] += rng.poisson(3.0, size=(epi.sum(), (gene_chrom == "chr7").sum()))
    counts[np.ix_(epi, gene_chrom == "chr10")] = rng.poisson(0.3, size=(epi.sum(), (gene_chrom == "chr10").sum()))

    adata = ad.AnnData(counts)
    adata.var_names = genes
    adata.var["chromosome"] = gene_chrom
    adata.var["start"] = [i % n_genes_per_chrom * 1000 + 1 for i in range(len(genes))]
    adata.var["end"] = adata.var["start"] + 500
    adata.obs_names = [f"cell{i}" for i in range(n_obs)]
    adata.obs["cell_id"] = adata.obs_names.to_numpy()
    adata.obs["cell_type"] = pd.Categorical(labels)
    adata.obsm["spatial"] = rng.random((n_obs, 2)) * 400.0
    adata.layers["raw_counts"] = adata.X.copy()
    return adata


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "synthetic_example.h5ad")
    build().write(out, compression="gzip")
    print(f"wrote {out}")
