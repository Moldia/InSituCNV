"""Build the small example dataset published as a GitHub release asset.

Usage::

    python scripts/make_example_dataset.py INPUT.h5ad -o insitucnv_example.h5ad \
        --reference-key cell_type --n-cells 4000

INPUT should be an image-based spatial AnnData with raw counts (in
``layers['raw_counts']`` or ``X``), ``obsm['spatial']`` and a cell-type column.
The script subsamples cells (stratified by the reference key), keeps the raw
counts / spatial / reference columns only, and writes a gzip-compressed h5ad plus
its sha256 so the value can be pasted into ``insitucnv/datasets.py``.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("insitucnv_example.h5ad"))
    parser.add_argument("--reference-key", default="cell_type")
    parser.add_argument("--n-cells", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    if "raw_counts" not in adata.layers:
        adata.layers["raw_counts"] = adata.X.copy()

    rng = np.random.default_rng(args.seed)
    if adata.n_obs > args.n_cells:
        groups = adata.obs[args.reference_key].astype(str)
        frac = args.n_cells / adata.n_obs
        keep = np.concatenate(
            [
                rng.choice(idx, size=max(1, int(round(len(idx) * frac))), replace=False)
                for idx in (np.where(groups == g)[0] for g in groups.unique())
            ]
        )
        adata = adata[np.sort(keep)].copy()

    slim = sc.AnnData(
        X=adata.layers["raw_counts"].copy(),
        obs=adata.obs[[args.reference_key]].copy(),
        var=adata.var[[]].copy(),
    )
    slim.layers["raw_counts"] = slim.X.copy()
    slim.obsm["spatial"] = np.asarray(adata.obsm["spatial"])
    slim.obs["cell_id"] = slim.obs_names.astype(str)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    slim.write(args.output, compression="gzip")

    digest = hashlib.sha256(args.output.read_bytes()).hexdigest()
    print(f"wrote {args.output} ({slim.n_obs} cells x {slim.n_vars} genes)")
    print(f"sha256: {digest}")


if __name__ == "__main__":
    main()
