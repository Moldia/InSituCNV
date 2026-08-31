from __future__ import annotations

import json

import scanpy as sc

from insitucnv import run_insitucnv
from insitucnv.tl.cnv import CNV_LEIDEN_PREFIX


def test_run_insitucnv_end_to_end(tiny_adata, gene_reference_csv, tmp_path):
    resolutions = [0.5, 1.0]
    result = run_insitucnv(
        tiny_adata,
        output_dir=tmp_path,
        reference_key="cell_type",
        reference_categories=["T_cells"],
        gene_reference_path=gene_reference_csv,
        smoothing_neighbors=15,
        window_size=10,
        step=5,
        cluster_resolutions=resolutions,
    )

    adata = result["adata"]
    assert "X_cnv" in adata.obsm
    assert result["cluster_keys"] == [f"{CNV_LEIDEN_PREFIX}{r:g}" for r in resolutions]
    assert result["primary_cluster_key"] in adata.obs

    assert (tmp_path / "adata_cnv.h5ad").exists()
    assert (tmp_path / "adata_input.h5ad").exists()
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    assert summary["reference_categories"] == ["T_cells"]
    assert summary["build_neighbors"] is True
    assert list((tmp_path / "plots").glob("*.png"))

    # round-trips
    sc.read_h5ad(tmp_path / "adata_cnv.h5ad")


def test_run_insitucnv_primary_resolution(tiny_adata, gene_reference_csv, tmp_path):
    result = run_insitucnv(
        tiny_adata,
        output_dir=tmp_path,
        reference_key="cell_type",
        reference_categories=["T_cells"],
        gene_reference_path=gene_reference_csv,
        smoothing_neighbors=15,
        window_size=10,
        step=5,
        cluster_resolutions=[0.5, 1.0],
        primary_resolution=1.0,
    )
    assert result["primary_cluster_key"] == f"{CNV_LEIDEN_PREFIX}1"
    assert result["primary_cluster_key"] in result["adata"].obs
