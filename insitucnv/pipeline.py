from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

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


def load_xenium_dataset(xenium_dir: str | Path, sample_id: str | None = None):
    import scanpy as sc

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
    import scanpy as sc

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
    import scanpy as sc

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


def run_xenium_cnv_protocol(
    adata,
    output_dir: str | Path,
    reference_key: str = "cell_type",
    reference_categories: list[str] | None = None,
    smoothing_neighbors: int = 100,
    window_size: int = 60,
    step: int = 10,
    lfc_clip: float = 4.0,
    cluster_resolutions: list[float] | None = None,
):
    """Run the Xenium CNV protocol on an annotated AnnData.

    Thin wrapper around :func:`insitucnv.workflow.run_insitucnv`: it evaluates the
    clustering-quality metrics and reports the metric-selected resolution as the
    primary CNV clustering. ``adata`` should already carry ``raw_counts``,
    ``spatial`` and a ``reference_key`` column (e.g. from
    :func:`load_xenium_dataset` + :func:`preprocess_expression` +
    :func:`annotate_cell_types`).

    Returns ``(adata, metrics, summary)``.
    """
    from insitucnv.workflow import run_insitucnv

    result = run_insitucnv(
        adata,
        output_dir=output_dir,
        reference_key=reference_key,
        reference_categories=reference_categories,
        smoothing_neighbors=smoothing_neighbors,
        window_size=window_size,
        step=step,
        lfc_clip=lfc_clip,
        cluster_resolutions=cluster_resolutions,
        evaluate_resolution_metrics=True,
        select_resolution_by_metrics=True,
    )
    out = result["adata"]
    summary = dict(result["summary"])
    if "sample" in out.obs.columns and out.n_obs:
        summary["sample_id"] = str(out.obs["sample"].iloc[0])
    return out, result["metrics"], summary
