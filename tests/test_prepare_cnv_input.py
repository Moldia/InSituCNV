from __future__ import annotations

import pandas as pd

from insitucnv.tl import prepare_cnv_input
from insitucnv.tl.cnv import _has_neighbor_graph


def test_prepare_cnv_input_builds_graph_when_absent(tiny_adata, gene_reference_csv):
    assert not _has_neighbor_graph(tiny_adata)

    out = prepare_cnv_input(
        tiny_adata,
        gene_reference_path=gene_reference_csv,
        smoothing_neighbors=15,
        copy=True,
    )

    assert _has_neighbor_graph(out)
    for layer in ("norm", "M", "log_norm"):
        assert layer in out.layers
    assert list(out.var.columns.intersection(["chromosome", "start", "end"])) == [
        "chromosome",
        "start",
        "end",
    ]
    assert out.var["chromosome"].notna().all()


def test_prepare_cnv_input_requires_graph_when_disabled(tiny_adata, gene_reference_csv):
    try:
        prepare_cnv_input(
            tiny_adata,
            gene_reference_path=gene_reference_csv,
            build_neighbors=False,
            copy=True,
        )
    except ValueError as exc:
        assert "build_neighbors=False" in str(exc)
    else:
        raise AssertionError("expected ValueError when graph is absent and build_neighbors=False")


def test_prepare_cnv_input_reuses_existing_graph(tiny_adata, gene_reference_csv):
    import scanpy as sc

    sc.pp.pca(tiny_adata, n_comps=10)
    sc.pp.neighbors(tiny_adata, n_neighbors=10)
    marker = tiny_adata.obsp["connectivities"].copy()

    out = prepare_cnv_input(tiny_adata, gene_reference_path=gene_reference_csv, smoothing_neighbors=10, copy=True)
    # graph was not rebuilt
    assert (out.obsp["connectivities"] != marker).nnz == 0
    assert isinstance(out.var, pd.DataFrame)
