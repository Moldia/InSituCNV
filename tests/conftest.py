from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

N_GENES_PER_CHROM = 20
CHROMS = ("chr1", "chr2")


@pytest.fixture
def gene_reference_csv(tmp_path_factory) -> str:
    """A local gene-coordinate table so tests never hit the network."""
    rows = []
    for chrom_idx, chrom in enumerate(CHROMS):
        for i in range(N_GENES_PER_CHROM):
            gene = f"g{chrom_idx * N_GENES_PER_CHROM + i}"
            start = i * 1000 + 1
            rows.append({"gene_name": gene, "chromosome": chrom, "start": start, "end": start + 500})
    path = tmp_path_factory.mktemp("ref") / "genes.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


@pytest.fixture
def tiny_adata() -> ad.AnnData:
    """~80 cells x 40 genes on 2 chromosomes; tumor cells carry a chr2 gain."""
    rng = np.random.default_rng(0)
    n_ref, n_tumor = 45, 45
    n_obs = n_ref + n_tumor
    n_genes = len(CHROMS) * N_GENES_PER_CHROM

    counts = rng.poisson(2.0, size=(n_obs, n_genes)).astype(np.float32)
    counts[n_ref:, N_GENES_PER_CHROM:] += rng.poisson(6.0, size=(n_tumor, N_GENES_PER_CHROM))

    adata = ad.AnnData(counts)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    adata.obs_names = [f"cell{i}" for i in range(n_obs)]
    adata.obs["cell_id"] = adata.obs_names.to_numpy()
    adata.obs["cell_type"] = pd.Categorical(["T_cells"] * n_ref + ["Epithelial"] * n_tumor)
    adata.obs["sample"] = "tiny"
    adata.obsm["spatial"] = rng.random((n_obs, 2)) * 100.0
    adata.layers["raw_counts"] = adata.X.copy()
    return adata
