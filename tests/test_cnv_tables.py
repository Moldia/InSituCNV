import anndata as ad
import numpy as np
import pandas as pd
import pytest

from insitucnv.tl import assign_cnv_status, assign_cnv_subclones, mean_cnv_per_gene


def _adata():
    obs = pd.DataFrame({"cluster": ["0", "1", "1"]}, index=["c1", "c2", "c3"])
    var = pd.DataFrame(
        {"chromosome": ["1", "1"], "start": [1, 2], "end": [10, 20]},
        index=["GENE1", "GENE2"],
    )
    data = ad.AnnData(X=np.ones((3, 2)), obs=obs, var=var)
    data.obsm["X_cnv"] = np.array([[0.0, 0.1], [1.0, 1.5], [1.2, 1.6]])
    data.layers["gene_values_cnv"] = data.obsm["X_cnv"].copy()
    return data


def test_assign_cnv_status_from_explicit_normal_clusters():
    data = _adata()
    assign_cnv_status(data, "cluster", normal_clusters=["0"])
    assert data.obs.loc["c1", "cnv_status"] == "normal"
    assert set(data.obs.loc[["c2", "c3"], "cnv_status"]) == {"tumor"}


def test_assign_cnv_status_requires_manual_cluster_selection():
    data = _adata()
    with pytest.raises(ValueError, match="No CNV status is inferred automatically"):
        assign_cnv_status(data, "cluster")


def test_assign_cnv_subclones_from_explicit_tumor_clusters():
    data = _adata()
    assign_cnv_subclones(data, "cluster", tumor_clusters=["1"], normal_clusters=["0"])
    assert data.obs.loc["c1", "cnv_subclones"] == "normal"
    assert set(data.obs.loc[["c2", "c3"], "cnv_subclones"]) == {"tumor_clone_1"}


def test_mean_cnv_per_gene_for_tumor_cells():
    data = _adata()
    data.obs["cnv_status"] = ["normal", "tumor", "tumor"]
    table = mean_cnv_per_gene(data)
    assert table["gene"].tolist() == ["GENE1", "GENE2"]
    np.testing.assert_allclose(table["CNV_value"], [1.1, 1.55])
