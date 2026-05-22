from __future__ import annotations

import json
from pathlib import Path

import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from insitucnv.analysis import find_optimal_clustering
from insitucnv.pp import add_genomic_positions
from insitucnv.tl import smooth_data_for_cnv

AUTO_MARKERS = {
    "T_cells": ["CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R"],
    "B_cells": ["CD79A", "MS4A1", "CD74", "HLA-DRA"],
    "Myeloid": ["LYZ", "C1QC", "FCER1G", "TYROBP"],
    "Plasma": ["JCHAIN", "MZB1", "SDC1", "XBP1"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN"],
    "Adipocytes": ["PLIN1", "FABP4", "ADIPOQ", "LEP"],
    "PVLs": ["RGS5", "MCAM", "CSPG4", "ACTA2", "PDGFRB"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17"],
}

REFERENCE_PRIORITY = ["T_cells", "B_cells", "Myeloid", "Plasma", "Fibroblast", "Endothelial", "Adipocytes", "PVLs"]


def load_xenium_dataset(xenium_dir: str | Path, sample_id: str | None = None):
    xenium_dir = Path(xenium_dir)
    matrix_path = xenium_dir / "cell_feature_matrix.h5"
    cells_path = xenium_dir / "cells.csv.gz"
    if not matrix_path.exists():
        raise FileNotFoundError(f"Could not find {matrix_path}.")
    if not cells_path.exists():
        raise FileNotFoundError(f"Could not find {cells_path}.")

    adata = sc.read_10x_h5(matrix_path)
    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.astype(str)

    cells = pd.read_csv(cells_path)
    if "cell_id" not in cells.columns:
        raise KeyError("Expected a 'cell_id' column in cells.csv.gz.")

    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    shared = adata.obs_names.intersection(cells.index)
    if shared.empty:
        raise ValueError("No overlapping cell IDs between cell_feature_matrix.h5 and cells.csv.gz.")

    adata = adata[shared].copy()
    adata.obs = adata.obs.join(cells.loc[shared], how="left")
    if {"x_centroid", "y_centroid"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    else:
        raise KeyError("cells.csv.gz must contain x_centroid and y_centroid columns.")

    adata.obs["cell_id"] = adata.obs_names.astype(str)
    adata.obs["sample"] = sample_id or xenium_dir.name
    adata.layers["raw_counts"] = adata.X.copy()
    return adata


def preprocess_expression(
    adata,
    min_counts: int = 20,
    min_genes: int = 10,
    n_neighbors: int = 15,
    leiden_resolution: float = 0.5,
):
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=5)
    if "raw_counts" not in adata.layers:
        adata.layers["raw_counts"] = adata.X.copy()

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars), subset=False)
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added="leiden")
    return adata


def _annotation_column_name(df: pd.DataFrame) -> str:
    for candidate in ("cell_type", "annotation"):
        if candidate in df.columns:
            return candidate
    raise KeyError("Annotation CSV must contain one of: cell_type, annotation.")


def _cell_id_column_name(df: pd.DataFrame) -> str:
    for candidate in ("cell_id", "cell.id"):
        if candidate in df.columns:
            return candidate
    raise KeyError("Annotation CSV must contain one of: cell_id, cell.id.")


def annotate_cell_types(
    adata,
    annotation_csv: str | Path | None = None,
    cluster_key: str = "leiden",
    min_score: float = 0.15,
    majority_threshold: float = 0.45,
    output_key: str = "cell_type",
):
    adata = adata.copy()

    if annotation_csv is not None:
        annotations = pd.read_csv(annotation_csv)
        cell_id_col = _cell_id_column_name(annotations)
        annotation_col = _annotation_column_name(annotations)
        annotations[cell_id_col] = annotations[cell_id_col].astype(str)
        lookup = annotations.set_index(cell_id_col)[annotation_col]
        adata.obs[output_key] = adata.obs["cell_id"].map(lookup).fillna("Unknown")
        return adata

    present_markers = {}
    for label, genes in AUTO_MARKERS.items():
        available = [gene for gene in genes if gene in adata.var_names]
        if len(available) >= 2:
            present_markers[label] = available

    if not present_markers:
        raise ValueError("None of the marker panels were found in the Xenium gene panel.")

    score_columns = []
    for label, genes in present_markers.items():
        score_name = f"{label}_score"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, use_raw=False)
        score_columns.append(score_name)

    score_df = adata.obs[score_columns].copy()
    adata.obs[f"{output_key}_score"] = score_df.max(axis=1)
    adata.obs[f"{output_key}_raw"] = score_df.idxmax(axis=1).str.replace("_score", "", regex=False)
    adata.obs[output_key] = np.where(
        adata.obs[f"{output_key}_score"] >= min_score,
        adata.obs[f"{output_key}_raw"],
        "Unknown",
    )

    if cluster_key in adata.obs:
        cluster_assignments = {}
        for cluster, sub_df in adata.obs.groupby(cluster_key):
            dominant = sub_df[output_key].value_counts(normalize=True)
            if dominant.empty:
                cluster_assignments[cluster] = "Unknown"
                continue
            label = dominant.index[0]
            frac = dominant.iloc[0]
            cluster_assignments[cluster] = label if frac >= majority_threshold else "Unknown"
        adata.obs[f"{output_key}_cluster_consensus"] = adata.obs[cluster_key].map(cluster_assignments).fillna("Unknown")
        adata.obs[output_key] = adata.obs[f"{output_key}_cluster_consensus"]

    return adata


def _resolve_reference_categories(adata, reference_key: str) -> list[str]:
    present = set(adata.obs[reference_key].astype(str).unique())
    selected = [label for label in REFERENCE_PRIORITY if label in present]
    if selected:
        return selected
    fallback = sorted(label for label in present if label not in {"Epithelial", "Unknown"})
    if fallback:
        return fallback
    raise ValueError(f"Could not determine inferCNV reference categories from {reference_key}.")


def _select_best_resolution(metrics: pd.DataFrame) -> float:
    if metrics.empty:
        raise ValueError("No valid clustering resolutions were evaluated.")

    ranked = metrics.copy()
    ranked["db_inverted"] = -ranked["davies_bouldin_score"]
    metric_cols = ["silhouette_score", "stability_score", "spatial_cohesion_score", "db_inverted"]
    for col in metric_cols:
        std = ranked[col].std()
        if pd.isna(std) or std == 0:
            ranked[f"{col}_z"] = 0.0
        else:
            ranked[f"{col}_z"] = (ranked[col] - ranked[col].mean()) / std
    ranked["combined_score"] = ranked[[f"{col}_z" for col in metric_cols]].sum(axis=1)
    best = ranked.sort_values(["combined_score", "stability_score", "silhouette_score"], ascending=False).iloc[0]
    return float(best["resolution"])


def _save_spatial_plot(adata, color: str, output_path: Path, title: str, size: float = 4.0):
    coords = adata.obsm["spatial"]
    values = adata.obs[color]

    fig, ax = plt.subplots(figsize=(8, 8))
    if pd.api.types.is_numeric_dtype(values):
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=size, cmap="viridis", linewidths=0)
        fig.colorbar(scatter, ax=ax, shrink=0.8)
    else:
        cats = pd.Categorical(values.astype(str))
        palette = plt.get_cmap("tab20", max(len(cats.categories), 1))
        for idx, category in enumerate(cats.categories):
            mask = cats == category
            ax.scatter(coords[mask, 0], coords[mask, 1], s=size, label=category, color=palette(idx), linewidths=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=3)

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.invert_yaxis()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def run_xenium_cnv_protocol(
    adata,
    output_dir: str | Path,
    reference_key: str = "cell_type",
    reference_categories: list[str] | None = None,
    smoothing_neighbors: int = 20,
    window_size: int = 60,
    step: int = 10,
    lfc_clip: float = 4.0,
    cluster_resolutions: list[float] | None = None,
):
    adata = adata.copy()
    output_dir = Path(output_dir)
    plots_dir = output_dir / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    adata.write(output_dir / "adata_preprocessed.h5ad", compression="gzip")

    adata.X = adata.layers["raw_counts"].copy()
    sc.pp.normalize_total(adata)
    adata.layers["norm"] = adata.X.copy()
    smooth_data_for_cnv(adata, layer_w_norm_counts="norm", n_neighbors=smoothing_neighbors)
    adata.layers["smooth_norm"] = adata.layers["M"].copy()

    adata.X = adata.layers["M"].copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata = add_genomic_positions(adata)

    reference_categories = reference_categories or _resolve_reference_categories(adata, reference_key)

    cnv.tl.infercnv(
        adata,
        window_size=window_size,
        step=step,
        lfc_clip=lfc_clip,
        reference_key=reference_key,
        reference_cat=reference_categories,
        chunksize=1000,
        calculate_gene_values=True,
    )
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)

    cluster_resolutions = cluster_resolutions or [0.1, 0.2, 0.3]
    metrics = find_optimal_clustering(adata, resolutions=cluster_resolutions)
    metrics.to_csv(output_dir / "cluster_resolution_metrics.csv", index=False)
    best_resolution = _select_best_resolution(metrics) if not metrics.empty else float(cluster_resolutions[0])
    cluster_key = f"cnv_leiden_{best_resolution:g}"
    cnv.tl.leiden(adata, resolution=best_resolution, key_added=cluster_key)

    _save_spatial_plot(adata, reference_key, plots_dir / "spatial_cell_types.png", "Spatial cell types")
    _save_spatial_plot(adata, cluster_key, plots_dir / "spatial_cnv_clusters.png", "Spatial CNV clusters")

    cnv.pl.chromosome_heatmap(adata, groupby=cluster_key, dendrogram=True, show=False, vmin=-0.4, vmax=0.4)
    plt.savefig(plots_dir / "cnv_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close()

    adata.write(output_dir / "adata_cnv.h5ad", compression="gzip")
    summary = {
        "sample_id": str(adata.obs["sample"].iloc[0]) if "sample" in adata.obs.columns and adata.n_obs else "unknown",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "reference_key": reference_key,
        "reference_categories": reference_categories,
        "best_resolution": best_resolution,
        "cluster_key": cluster_key,
        "smoothing_neighbors": smoothing_neighbors,
        "window_size": window_size,
        "step": step,
        "lfc_clip": lfc_clip,
    }
    (output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    return adata, metrics, summary
