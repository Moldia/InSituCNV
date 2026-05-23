from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd
from scipy import sparse


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
    import scanpy as sc

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
    import scanpy as sc

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
    from insitucnv.tl.moments import smooth_data_for_cnv

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
        from insitucnv.pp import add_genomic_positions

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
    import infercnvpy as cnv

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
    copy: bool = False,
    **neighbors_kwargs: Any,
):
    """Compute the CNV PCA and neighbor graph used for clustering."""
    import infercnvpy as cnv

    out = adata.copy() if copy else adata
    if run_pca:
        cnv.tl.pca(out)
    cnv.pp.neighbors(out, **neighbors_kwargs)
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
    from matplotlib import cm
    from matplotlib.colors import to_hex

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
    import infercnvpy as cnv

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
                import scanpy as sc

                sc.tl.dendrogram(out, groupby=key)
            except Exception as exc:  # pragma: no cover - plotting convenience should not stop analysis
                warnings.warn(f"Could not compute dendrogram for {key}: {exc}", RuntimeWarning, stacklevel=2)
        keys.append(key)

    out.uns.setdefault("insitucnv", {})["cluster_keys"] = keys
    return keys


def assign_cnv_status(
    adata,
    cluster_key: str,
    tumor_clusters: Sequence[str] | None = None,
    normal_clusters: Sequence[str] | None = None,
    output_key: str = "cnv_status",
    tumor_label: str = "tumor",
    normal_label: str = "normal",
    unassigned_label: str = "unassigned",
):
    """Annotate CNV clusters as tumor/normal from explicit cluster choices."""
    tumor_values = [str(cluster) for cluster in tumor_clusters or []]
    normal_values = [str(cluster) for cluster in normal_clusters or []]
    if not tumor_values and not normal_values:
        raise ValueError("Pass tumor_clusters, normal_clusters, or both. No CNV status is inferred automatically.")

    labels = adata.obs[cluster_key].astype(str)
    all_clusters = set(labels.unique())
    tumor_set = set(tumor_values)
    normal_set = set(normal_values)

    overlap = tumor_set & normal_set
    if overlap:
        raise ValueError(f"Clusters cannot be both tumor and normal: {sorted(overlap)}")
    unknown = (tumor_set | normal_set) - all_clusters
    if unknown:
        raise ValueError(f"Clusters not found in adata.obs['{cluster_key}']: {sorted(unknown)}")

    if not tumor_values:
        tumor_set = all_clusters - normal_set
    if not normal_values:
        normal_set = all_clusters - tumor_set

    status = pd.Series(unassigned_label, index=adata.obs_names, dtype="object")
    status.loc[labels.isin(normal_set).to_numpy()] = normal_label
    status.loc[labels.isin(tumor_set).to_numpy()] = tumor_label

    categories = [normal_label, tumor_label]
    if unassigned_label in set(status):
        categories.append(unassigned_label)
    adata.obs[output_key] = pd.Categorical(status, categories=categories)
    return adata


def assign_cnv_subclones(
    adata,
    cluster_key: str,
    tumor_clusters: Sequence[str],
    output_key: str = "cnv_subclones",
    normal_clusters: Sequence[str] | None = None,
    normal_label: str = "normal",
    other_tumor_label: str = "tumor",
    prefix: str = "tumor_clone_",
):
    """Label explicitly selected tumor CNV clusters as separate clone groups."""
    labels = adata.obs[cluster_key].astype(str)
    all_clusters = set(labels.unique())
    tumor_set = {str(cluster) for cluster in tumor_clusters}
    normal_set = {str(cluster) for cluster in normal_clusters or []}

    overlap = tumor_set & normal_set
    if overlap:
        raise ValueError(f"Clusters cannot be both tumor clones and normal: {sorted(overlap)}")
    unknown = (tumor_set | normal_set) - all_clusters
    if unknown:
        raise ValueError(f"Clusters not found in adata.obs['{cluster_key}']: {sorted(unknown)}")

    if normal_clusters is None:
        normal_set = all_clusters - tumor_set

    def _label(cluster):
        cluster = str(cluster)
        if cluster in tumor_set:
            return f"{prefix}{cluster}"
        if cluster in normal_set:
            return normal_label
        return other_tumor_label

    values = labels.map(_label)
    categories = [normal_label, *[f"{prefix}{cluster}" for cluster in sorted(tumor_set)]]
    if other_tumor_label in set(values):
        categories.append(other_tumor_label)
    adata.obs[output_key] = pd.Categorical(values, categories=categories)
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
