from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any
import warnings

import infercnvpy as cnv
from matplotlib import cm
from matplotlib.colors import to_hex
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from insitucnv.pp import add_genomic_positions
from insitucnv.tl.moments import smooth_data_for_cnv


def _copy_matrix(matrix):
    return matrix.copy() if hasattr(matrix, "copy") else matrix


def _set_x_from_layer(adata, layer: str | None):
    if layer is None:
        return
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")
    adata.X = _copy_matrix(adata.layers[layer])


def _as_dense_vector(values) -> np.ndarray:
    return np.asarray(values).ravel()


def _resolution_label(resolution: float) -> str:
    return f"{float(resolution):g}"


def _cluster_key(prefix: str, resolution: float) -> str:
    return f"{prefix}{_resolution_label(resolution)}"


def normalize_counts(
    adata,
    input_layer: str | None = "raw_counts",
    output_layer: str = "norm",
    target_sum: float | None = 1e4,
    copy: bool = False,
):
    """Normalize counts and store them in a layer.

    ``adata.X`` is set to the normalized matrix, matching the notebook workflow.
    The modified AnnData object is always returned.
    """
    out = adata.copy() if copy else adata
    _set_x_from_layer(out, input_layer)
    sc.pp.normalize_total(out, target_sum=target_sum)
    out.layers[output_layer] = _copy_matrix(out.X)
    out.uns.setdefault("insitucnv", {})["normalization"] = {
        "input_layer": input_layer,
        "output_layer": output_layer,
        "target_sum": target_sum,
    }
    return out


def log_normalize_counts(
    adata,
    input_layer: str = "M",
    output_layer: str = "log_norm",
    target_sum: float | None = 1e4,
    copy: bool = False,
):
    """Normalize and log-transform a layer for ``infercnvpy`` input."""
    out = adata.copy() if copy else adata
    _set_x_from_layer(out, input_layer)
    sc.pp.normalize_total(out, target_sum=target_sum)
    sc.pp.log1p(out)
    out.layers[output_layer] = _copy_matrix(out.X)
    out.uns.setdefault("insitucnv", {})["log_normalization"] = {
        "input_layer": input_layer,
        "output_layer": output_layer,
        "target_sum": target_sum,
    }
    return out


def prepare_cnv_input(
    adata,
    raw_layer: str = "raw_counts",
    normalized_layer: str = "norm",
    smoothed_layer: str = "M",
    log_layer: str = "log_norm",
    target_sum: float | None = 1e4,
    smoothing_neighbors: int = 100,
    smoothing_mode: str = "connectivities",
    add_gene_positions: bool = True,
    gene_reference: pd.DataFrame | None = None,
    gene_reference_path: str | Path | None = None,
    drop_unmapped_genes: bool = True,
    copy: bool = False,
):
    """Run the standard InSituCNV preprocessing chain.

    Steps:
    1. restore raw counts from ``raw_layer``;
    2. normalize counts into ``normalized_layer``;
    3. smooth normalized counts over the neighbor graph into ``smoothed_layer``;
    4. normalize/log-transform smoothed counts into ``log_layer`` and ``adata.X``;
    5. optionally add genomic positions required by ``infercnvpy``.
    """
    out = adata.copy() if copy else adata
    if raw_layer not in out.layers:
        out.layers[raw_layer] = _copy_matrix(out.X)

    normalize_counts(out, input_layer=raw_layer, output_layer=normalized_layer, target_sum=target_sum)
    smooth_data_for_cnv(
        out,
        layer_w_norm_counts=normalized_layer,
        n_neighbors=smoothing_neighbors,
        mode=smoothing_mode,
        output_layer=smoothed_layer,
    )
    log_normalize_counts(out, input_layer=smoothed_layer, output_layer=log_layer, target_sum=target_sum)

    if add_gene_positions:
        out = add_genomic_positions(
            out,
            reference=gene_reference,
            reference_path=gene_reference_path,
            drop_unmapped_genes=drop_unmapped_genes,
        )

    out.uns.setdefault("insitucnv", {})["cnv_input"] = {
        "raw_layer": raw_layer,
        "normalized_layer": normalized_layer,
        "smoothed_layer": smoothed_layer,
        "log_layer": log_layer,
        "target_sum": target_sum,
        "smoothing_neighbors": smoothing_neighbors,
        "smoothing_mode": smoothing_mode,
    }
    return out


def run_infercnv(
    adata,
    reference_key: str,
    reference_categories: Sequence[str] | None = None,
    input_layer: str | None = "log_norm",
    window_size: int = 60,
    step: int = 10,
    lfc_clip: float = 4.0,
    chunksize: int = 1000,
    calculate_gene_values: bool = True,
    copy: bool = False,
    **kwargs: Any,
):
    """Run ``infercnvpy.tl.infercnv`` with manuscript-style defaults."""
    out = adata.copy() if copy else adata
    _set_x_from_layer(out, input_layer)
    infercnv_kwargs = {
        "window_size": window_size,
        "step": step,
        "lfc_clip": lfc_clip,
        "reference_key": reference_key,
        "chunksize": chunksize,
        "calculate_gene_values": calculate_gene_values,
        **kwargs,
    }
    if reference_categories is not None:
        infercnv_kwargs["reference_cat"] = list(reference_categories)
    cnv.tl.infercnv(out, **infercnv_kwargs)
    out.uns.setdefault("insitucnv", {})["infercnv"] = {
        "reference_key": reference_key,
        "reference_categories": list(reference_categories) if reference_categories is not None else None,
        "input_layer": input_layer,
        "window_size": window_size,
        "step": step,
        "lfc_clip": lfc_clip,
        "chunksize": chunksize,
        "calculate_gene_values": calculate_gene_values,
    }
    return out


def compute_cnv_neighbors(
    adata,
    run_pca: bool = True,
    run_umap: bool = False,
    copy: bool = False,
    **neighbors_kwargs: Any,
):
    """Compute the CNV PCA and neighbor graph used for clustering."""
    out = adata.copy() if copy else adata
    if run_pca:
        cnv.tl.pca(out)
    cnv.pp.neighbors(out, **neighbors_kwargs)
    if run_umap:
        cnv.tl.umap(out)
    return out


def _subset_mask(adata, subset_key: str | None, subset_values: Sequence[str] | None):
    if subset_key is None:
        return None
    if subset_key not in adata.obs:
        raise KeyError(f"obs column '{subset_key}' not found.")
    if subset_values is None:
        raise ValueError("subset_values must be provided when subset_key is used.")
    return adata.obs[subset_key].astype(str).isin({str(value) for value in subset_values}).to_numpy()


def _set_cluster_palette(adata, key: str, outside_label: str | None = None, palette: str = "tab20"):
    cats = pd.Categorical(adata.obs[key].astype(str)).categories
    cmap = cm.get_cmap(palette, max(len(cats), 1))
    colors = []
    for idx, category in enumerate(cats):
        if outside_label is not None and str(category) == str(outside_label):
            colors.append("#D3D3D3")
        else:
            colors.append(to_hex(cmap(idx)))
    adata.obs[key] = pd.Categorical(adata.obs[key].astype(str), categories=cats)
    adata.uns[f"{key}_colors"] = colors


def cluster_cnv_resolutions(
    adata,
    resolutions: Sequence[float],
    key_prefix: str = "cnv_leiden_res",
    subset_key: str | None = None,
    subset_values: Sequence[str] | None = None,
    subset_label_prefix: str = "",
    outside_label: str = "other",
    dendrogram: bool = False,
    palette: str = "tab20",
    copy: bool = False,
    **leiden_kwargs: Any,
) -> list[str]:
    """Cluster CNV profiles at one or more Leiden resolutions.

    For epithelial-only or tumor-only clustering, pass ``subset_key`` and
    ``subset_values``. The subset labels are written back to the full AnnData
    object, with non-subset cells assigned ``outside_label``.
    """
    out = adata.copy() if copy else adata
    mask = _subset_mask(out, subset_key, subset_values)
    keys = []

    for resolution in resolutions:
        key = _cluster_key(key_prefix, resolution)
        if mask is None:
            cnv.tl.leiden(out, resolution=float(resolution), key_added=key, **leiden_kwargs)
            outside_for_palette = None
        else:
            subset = out[mask].copy()
            cnv.tl.leiden(subset, resolution=float(resolution), key_added=key, **leiden_kwargs)
            out.obs[key] = outside_label
            out.obs.loc[mask, key] = subset_label_prefix + subset.obs[key].astype(str).to_numpy()
            outside_for_palette = outside_label

        _set_cluster_palette(out, key, outside_label=outside_for_palette, palette=palette)
        if dendrogram:
            try:
                sc.tl.dendrogram(out, groupby=key)
            except Exception as exc:  # pragma: no cover - plotting convenience should not stop analysis
                warnings.warn(f"Could not compute dendrogram for {key}: {exc}", RuntimeWarning, stacklevel=2)
        keys.append(key)

    out.uns.setdefault("insitucnv", {})["cluster_keys"] = keys
    return keys


def calculate_cnv_burden(
    adata,
    matrix_key: str = "X_cnv",
    output_key: str = "cnv_burden",
) -> pd.Series:
    """Calculate per-cell mean absolute CNV signal."""
    if matrix_key in adata.obsm:
        matrix = adata.obsm[matrix_key]
    elif matrix_key in adata.layers:
        matrix = adata.layers[matrix_key]
    else:
        raise KeyError(f"Could not find '{matrix_key}' in adata.obsm or adata.layers.")

    if sparse.issparse(matrix):
        burden = _as_dense_vector(abs(matrix).mean(axis=1))
    else:
        burden = np.mean(np.abs(np.asarray(matrix)), axis=1)
    adata.obs[output_key] = burden
    return adata.obs[output_key]


def assign_cnv_status(
    adata,
    cluster_key: str,
    tumor_clusters: Sequence[str] | None = None,
    normal_clusters: Sequence[str] | None = None,
    output_key: str = "cnv_status",
    tumor_label: str = "tumor",
    normal_label: str = "normal",
    burden_key: str = "cnv_burden",
    matrix_key: str = "X_cnv",
):
    """Annotate CNV clusters as tumor/normal.

    If neither ``tumor_clusters`` nor ``normal_clusters`` is supplied, the
    cluster with the lowest CNV burden is treated as normal-like.
    """
    labels = adata.obs[cluster_key].astype(str)
    if tumor_clusters is None and normal_clusters is None:
        if burden_key not in adata.obs:
            calculate_cnv_burden(adata, matrix_key=matrix_key, output_key=burden_key)
        normal_cluster = adata.obs.groupby(cluster_key, observed=False)[burden_key].mean().idxmin()
        normal_clusters = [str(normal_cluster)]

    if tumor_clusters is not None:
        tumor_set = {str(cluster) for cluster in tumor_clusters}
        adata.obs[output_key] = np.where(labels.isin(tumor_set), tumor_label, normal_label)
    else:
        normal_set = {str(cluster) for cluster in normal_clusters or []}
        adata.obs[output_key] = np.where(labels.isin(normal_set), normal_label, tumor_label)

    adata.obs[output_key] = pd.Categorical(adata.obs[output_key], categories=[tumor_label, normal_label])
    return adata


def assign_cnv_subclones(
    adata,
    cluster_key: str,
    tumor_clusters: Sequence[str],
    output_key: str = "cnv_subclones",
    normal_label: str = "normal",
    prefix: str = "subclone_",
):
    """Label selected tumor CNV clusters as subclones and all others as normal."""
    tumor_set = {str(cluster) for cluster in tumor_clusters}

    def _label(cluster):
        cluster = str(cluster)
        return f"{prefix}{cluster}" if cluster in tumor_set else normal_label

    adata.obs[output_key] = pd.Categorical(adata.obs[cluster_key].astype(str).map(_label))
    return adata


def mean_cnv_per_gene(
    adata,
    layer: str = "gene_values_cnv",
    obs_key: str | None = "cnv_status",
    obs_values: Sequence[str] | None = ("tumor",),
    mask: Sequence[bool] | None = None,
) -> pd.DataFrame:
    """Return mean CNV value per gene for a selected set of cells."""
    if layer is None:
        matrix = adata.X
    elif layer in adata.layers:
        matrix = adata.layers[layer]
    else:
        raise KeyError(f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")

    if mask is None:
        if obs_key is None:
            mask_arr = np.ones(adata.n_obs, dtype=bool)
        else:
            if obs_key not in adata.obs:
                raise KeyError(f"obs column '{obs_key}' not found.")
            values = {str(value) for value in (obs_values or [])}
            mask_arr = adata.obs[obs_key].astype(str).isin(values).to_numpy()
    else:
        mask_arr = np.asarray(mask, dtype=bool)

    if mask_arr.shape[0] != adata.n_obs:
        raise ValueError("mask length must match adata.n_obs.")
    if not mask_arr.any():
        raise ValueError("No cells selected for mean CNV calculation.")

    selected = matrix[mask_arr, :]
    if sparse.issparse(selected):
        mean_values = _as_dense_vector(selected.mean(axis=0))
    else:
        mean_values = np.asarray(selected).mean(axis=0)

    result = pd.DataFrame({"gene": adata.var_names.astype(str), "CNV_value": mean_values})
    for col in ("chromosome", "start", "end"):
        if col in adata.var:
            result[col] = adata.var[col].to_numpy()
    ordered = ["gene", *[col for col in ("chromosome", "start", "end") if col in result], "CNV_value"]
    return result.loc[:, ordered]


def export_mean_cnv_per_gene(
    adata,
    output_path: str | Path,
    layer: str = "gene_values_cnv",
    obs_key: str | None = "cnv_status",
    obs_values: Sequence[str] | None = ("tumor",),
    sep: str = "\t",
) -> pd.DataFrame:
    """Write ``mean_cnv_per_gene`` output and return the table."""
    table = mean_cnv_per_gene(adata, layer=layer, obs_key=obs_key, obs_values=obs_values)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, sep=sep, index=False)
    return table


def export_cell_groups(
    adata,
    output_path: str | Path,
    group_key: str,
    cell_id_key: str = "cell_id",
    group_col: str = "group",
) -> pd.DataFrame:
    """Export a Xenium Explorer-compatible cell group CSV."""
    if group_key not in adata.obs:
        raise KeyError(f"obs column '{group_key}' not found.")
    if cell_id_key in adata.obs:
        cell_ids = adata.obs[cell_id_key].astype(str)
    else:
        cell_ids = adata.obs_names.astype(str)

    table = pd.DataFrame({"cell_id": cell_ids.to_numpy(), group_col: adata.obs[group_key].astype(str).to_numpy()})
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    return table
