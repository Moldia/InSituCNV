from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from insitucnv.tl import prepare_cnv_input, smooth_data_for_cnv
from insitucnv.tl.cnv import SMOOTHING_NEIGHBORS_KEY, _graph_n_neighbors


def _expression_graph(adata, n_neighbors):
    import scanpy as sc

    sc.pp.pca(adata, n_comps=10)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)


def test_prepare_cnv_input_builds_graph_when_absent(tiny_adata, gene_reference_csv):
    assert _graph_n_neighbors(tiny_adata, "neighbors") is None

    out = prepare_cnv_input(
        tiny_adata,
        gene_reference_path=gene_reference_csv,
        smoothing_neighbors=15,
        copy=True,
    )

    # the smoothing graph is built with exactly smoothing_neighbors
    assert _graph_n_neighbors(out, SMOOTHING_NEIGHBORS_KEY) == 15
    assert out.uns["insitucnv"]["cnv_input"]["smoothing_neighbors_key"] == SMOOTHING_NEIGHBORS_KEY
    for layer in ("norm", "M", "log_norm"):
        assert layer in out.layers
    assert list(out.var.columns.intersection(["chromosome", "start", "end"])) == [
        "chromosome",
        "start",
        "end",
    ]
    assert out.var["chromosome"].notna().all()


def test_prepare_cnv_input_requires_graph_when_disabled(tiny_adata, gene_reference_csv):
    with pytest.raises(ValueError, match="build_neighbors=False"):
        prepare_cnv_input(
            tiny_adata,
            gene_reference_path=gene_reference_csv,
            build_neighbors=False,
            copy=True,
        )


def test_prepare_cnv_input_rejects_mismatched_graph_when_disabled(tiny_adata, gene_reference_csv):
    _expression_graph(tiny_adata, n_neighbors=10)
    with pytest.raises(ValueError, match="n_neighbors=10"):
        prepare_cnv_input(
            tiny_adata,
            gene_reference_path=gene_reference_csv,
            smoothing_neighbors=15,
            build_neighbors=False,
            copy=True,
        )


def test_prepare_cnv_input_reuses_matching_graph(tiny_adata, gene_reference_csv):
    _expression_graph(tiny_adata, n_neighbors=10)
    marker = tiny_adata.obsp["connectivities"].copy()

    out = prepare_cnv_input(tiny_adata, gene_reference_path=gene_reference_csv, smoothing_neighbors=10, copy=True)
    # graph was not rebuilt
    assert (out.obsp["connectivities"] != marker).nnz == 0
    assert SMOOTHING_NEIGHBORS_KEY not in out.uns
    assert out.uns["insitucnv"]["cnv_input"]["smoothing_neighbors_key"] == "neighbors"
    assert isinstance(out.var, pd.DataFrame)


def test_prepare_cnv_input_builds_own_graph_on_mismatch(tiny_adata, gene_reference_csv):
    _expression_graph(tiny_adata, n_neighbors=10)
    marker_graph = tiny_adata.obsp["connectivities"].copy()
    marker_pca = tiny_adata.obsm["X_pca"].copy()

    out = prepare_cnv_input(tiny_adata, gene_reference_path=gene_reference_csv, smoothing_neighbors=15, copy=True)

    assert _graph_n_neighbors(out, SMOOTHING_NEIGHBORS_KEY) == 15
    # the user's expression graph and PCA are left untouched
    assert _graph_n_neighbors(out, "neighbors") == 10
    assert (out.obsp["connectivities"] != marker_graph).nnz == 0
    np.testing.assert_array_equal(out.obsm["X_pca"], marker_pca)


def test_smooth_data_for_cnv_rejects_mismatched_graph(tiny_adata):
    _expression_graph(tiny_adata, n_neighbors=10)
    with pytest.raises(ValueError, match="n_neighbors=10"):
        smooth_data_for_cnv(tiny_adata, n_neighbors=20)


def test_smooth_data_for_cnv_matches_scvelo(tiny_adata):
    """With a graph built with k, smoothing is unchanged from the previous scvelo-based code."""
    from scipy.sparse import csr_matrix
    from scvelo.preprocessing.neighbors import get_connectivities

    _expression_graph(tiny_adata, n_neighbors=10)
    smooth_data_for_cnv(tiny_adata, n_neighbors=10)

    weights = get_connectivities(tiny_adata, "connectivities", n_neighbors=10, recurse_neighbors=False)
    expected = csr_matrix.dot(weights, csr_matrix(tiny_adata.X)).toarray()
    np.testing.assert_allclose(tiny_adata.layers["M"], expected, rtol=1e-5)
