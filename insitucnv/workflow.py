from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from insitucnv.pl import plot_chromosome_heatmap, plot_spatial
from insitucnv.tl import (
    cluster_cnv_resolutions,
    compute_cnv_neighbors,
    prepare_cnv_input,
    run_infercnv,
)
from insitucnv.tl.cnv import cnv_leiden_key

REFERENCE_PRIORITY = ["T_cells", "B_cells", "Myeloid", "Plasma", "Fibroblast", "Endothelial", "Adipocytes", "PVLs"]


def resolve_reference_categories(adata, reference_key: str, priority: list[str] | None = None) -> list[str]:
    """Choose non-tumor reference categories for inferCNV."""
    if reference_key not in adata.obs:
        raise KeyError(f"adata.obs['{reference_key}'] not found.")
    present = set(adata.obs[reference_key].astype(str).unique())
    selected = [label for label in (priority or REFERENCE_PRIORITY) if label in present]
    if selected:
        return selected
    fallback = sorted(label for label in present if label.lower() not in {"epithelial", "tumor", "unknown"})
    if fallback:
        warnings.warn(
            "No priority reference cell types were found in "
            f"adata.obs['{reference_key}']. Falling back to every category except "
            f"epithelial/tumor/unknown as the inferCNV reference: {fallback}. "
            "Pass reference_categories explicitly to control this.",
            RuntimeWarning,
            stacklevel=2,
        )
        return fallback
    raise ValueError(f"Could not determine inferCNV reference categories from '{reference_key}'.")


def select_best_resolution(metrics: pd.DataFrame) -> float:
    """Select a Leiden resolution from the metrics returned by ``find_optimal_clustering``."""
    if metrics.empty:
        raise ValueError("No valid clustering resolutions were evaluated.")

    ranked = metrics.copy()
    ranked["db_inverted"] = -ranked["davies_bouldin_score"]
    metric_cols = ["silhouette_score", "stability_score", "spatial_cohesion_score", "db_inverted"]
    for col in metric_cols:
        std = ranked[col].std()
        ranked[f"{col}_z"] = 0.0 if pd.isna(std) or std == 0 else (ranked[col] - ranked[col].mean()) / std
    ranked["combined_score"] = ranked[[f"{col}_z" for col in metric_cols]].sum(axis=1)
    best = ranked.sort_values(["combined_score", "stability_score", "silhouette_score"], ascending=False).iloc[0]
    return float(best["resolution"])


def _json_ready(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def run_insitucnv(
    adata,
    output_dir: str | Path,
    reference_key: str,
    reference_categories: list[str] | None = None,
    raw_layer: str = "raw_counts",
    target_sum: float | None = 1e4,
    smoothing_neighbors: int = 100,
    smoothing_mode: str = "connectivities",
    build_neighbors: bool = True,
    neighbors_n_neighbors: int = 15,
    neighbors_n_pcs: int = 50,
    gene_reference_path: str | Path | None = None,
    window_size: int = 60,
    step: int = 10,
    lfc_clip: float = 4.0,
    chunksize: int = 1000,
    cluster_resolutions: list[float] | None = None,
    primary_resolution: float | None = None,
    select_resolution_by_metrics: bool = False,
    evaluate_resolution_metrics: bool = False,
    spatial_key: str = "spatial",
    point_size: float = 4.0,
    save_intermediate: bool = True,
    copy: bool = True,
    **infercnv_kwargs: Any,
) -> dict[str, Any]:
    """Run the core InSituCNV workflow on a prepared AnnData object.

    The input AnnData should contain raw counts in ``raw_layer`` or in ``X``,
    reference labels in ``adata.obs[reference_key]``, and spatial coordinates in
    ``adata.obsm[spatial_key]`` if spatial plots are requested. A neighbor graph
    for smoothing is built automatically when absent (``build_neighbors=True``);
    pass ``build_neighbors=False`` to require a precomputed graph. Gene symbols in
    ``var_names`` are matched to genomic coordinates using the bundled
    ``infercnvpy`` human (GRCh38) table unless ``gene_reference_path`` points to a
    CSV/TSV with ``gene_name,chromosome,start,end`` columns.

    Tumor/normal and clone labels are intentionally left to downstream manual
    annotation after the heatmap and spatial plots are reviewed.
    """
    out = adata.copy() if copy else adata
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    cluster_resolutions = cluster_resolutions or [0.1, 0.2, 0.3]
    if save_intermediate:
        out.write(output_dir / "adata_input.h5ad", compression="gzip")

    out = prepare_cnv_input(
        out,
        raw_layer=raw_layer,
        target_sum=target_sum,
        smoothing_neighbors=smoothing_neighbors,
        smoothing_mode=smoothing_mode,
        build_neighbors=build_neighbors,
        neighbors_n_neighbors=neighbors_n_neighbors,
        neighbors_n_pcs=neighbors_n_pcs,
        gene_reference_path=gene_reference_path,
        copy=False,
    )
    reference_categories = reference_categories or resolve_reference_categories(out, reference_key)

    run_infercnv(
        out,
        reference_key=reference_key,
        reference_categories=reference_categories,
        window_size=window_size,
        step=step,
        lfc_clip=lfc_clip,
        chunksize=chunksize,
        **infercnv_kwargs,
    )
    compute_cnv_neighbors(out)
    cluster_keys = cluster_cnv_resolutions(out, cluster_resolutions, dendrogram=True)

    metrics = pd.DataFrame()
    if evaluate_resolution_metrics or select_resolution_by_metrics:
        from insitucnv.analysis import find_optimal_clustering

        metrics = find_optimal_clustering(out, resolutions=cluster_resolutions, spatial_key=spatial_key)
        metrics.to_csv(output_dir / "cluster_resolution_metrics.csv", index=False)

    if select_resolution_by_metrics and not metrics.empty:
        primary_resolution = select_best_resolution(metrics)
    if primary_resolution is not None:
        primary_cluster_key = cnv_leiden_key(primary_resolution)
    else:
        primary_cluster_key = cluster_keys[0]

    for key in cluster_keys:
        plot_chromosome_heatmap(out, groupby=key, output_path=plots_dir / f"{key}_heatmap.png")
        if spatial_key in out.obsm:
            plot_spatial(
                out,
                color=key,
                output_path=plots_dir / f"{key}_spatial.png",
                spatial_key=spatial_key,
                point_size=point_size,
                title=f"Spatial {key}",
            )

    output_h5ad = output_dir / "adata_cnv.h5ad"
    out.write(output_h5ad, compression="gzip")

    summary = {
        "n_cells": int(out.n_obs),
        "n_genes": int(out.n_vars),
        "reference_key": reference_key,
        "reference_categories": reference_categories,
        "raw_layer": raw_layer,
        "target_sum": target_sum,
        "smoothing_neighbors": smoothing_neighbors,
        "smoothing_mode": smoothing_mode,
        "build_neighbors": build_neighbors,
        "window_size": window_size,
        "step": step,
        "lfc_clip": lfc_clip,
        "cluster_resolutions": cluster_resolutions,
        "cluster_keys": cluster_keys,
        "primary_cluster_key": primary_cluster_key,
        "output_h5ad": output_h5ad,
    }
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2, default=_json_ready))

    return {
        "adata": out,
        "cluster_keys": cluster_keys,
        "primary_cluster_key": primary_cluster_key,
        "metrics": metrics,
        "summary": summary,
    }
