from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def _prepare_output_path(output_path: str | Path | None):
    if output_path is None:
        return None
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def plot_chromosome_heatmap(
    adata,
    groupby: str,
    output_path: str | Path | None = None,
    dendrogram: bool = True,
    vmin: float = -0.4,
    vmax: float = 0.4,
    show: bool = False,
    **kwargs,
):
    """Plot and optionally save an ``infercnvpy`` chromosome heatmap."""
    cnv.pl.chromosome_heatmap(
        adata,
        groupby=groupby,
        dendrogram=dendrogram,
        vmin=vmin,
        vmax=vmax,
        show=show,
        **kwargs,
    )
    output_path = _prepare_output_path(output_path)
    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close()


def plot_embedding(
    adata,
    color: str | Sequence[str],
    basis: str = "umap",
    output_path: str | Path | None = None,
    show: bool = False,
    **kwargs,
):
    """Plot a Scanpy embedding such as UMAP colored by CNV clusters."""
    sc.pl.embedding(adata, basis=basis, color=color, show=show, **kwargs)
    output_path = _prepare_output_path(output_path)
    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close()


def plot_spatial(
    adata,
    color: str,
    output_path: str | Path | None = None,
    spatial_key: str = "spatial",
    point_size: float = 4.0,
    palette: str = "tab20",
    invert_yaxis: bool = True,
    title: str | None = None,
    show: bool = False,
):
    """Plot spatial coordinates colored by an ``adata.obs`` column."""
    if spatial_key not in adata.obsm:
        raise KeyError(f"adata.obsm['{spatial_key}'] not found.")
    if color not in adata.obs:
        raise KeyError(f"adata.obs['{color}'] not found.")

    coords = np.asarray(adata.obsm[spatial_key])
    values = adata.obs[color]

    fig, ax = plt.subplots(figsize=(8, 8))
    if pd.api.types.is_numeric_dtype(values):
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=point_size, cmap="viridis", linewidths=0)
        fig.colorbar(scatter, ax=ax, shrink=0.8)
    else:
        categories = pd.Categorical(values.astype(str))
        cmap = plt.get_cmap(palette, max(len(categories.categories), 1))
        for idx, category in enumerate(categories.categories):
            mask = categories == category
            ax.scatter(coords[mask, 0], coords[mask, 1], s=point_size, color=cmap(idx), label=category, linewidths=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=3)

    ax.set_title(title or color)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if invert_yaxis:
        ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")

    output_path = _prepare_output_path(output_path)
    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_cluster_composition(
    adata,
    cluster_key: str,
    annotation_key: str,
    output_path: str | Path | None = None,
    normalize: str = "index",
    figsize: tuple[float, float] = (8, 5),
):
    """Plot annotation proportions for each CNV cluster."""
    if cluster_key not in adata.obs:
        raise KeyError(f"adata.obs['{cluster_key}'] not found.")
    if annotation_key not in adata.obs:
        raise KeyError(f"adata.obs['{annotation_key}'] not found.")

    table = pd.crosstab(adata.obs[cluster_key], adata.obs[annotation_key], normalize=normalize)
    colors = None
    color_key = f"{annotation_key}_colors"
    if color_key in adata.uns and hasattr(adata.obs[annotation_key], "cat"):
        categories = list(adata.obs[annotation_key].cat.categories)
        color_lookup = dict(zip(categories, adata.uns[color_key]))
        colors = [color_lookup[col] for col in table.columns if col in color_lookup]

    ax = table.plot(kind="bar", stacked=True, figsize=figsize, color=colors, width=0.8)
    ax.set_ylabel("Proportion" if normalize else "Count")
    ax.set_xlabel("CNV cluster")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title=annotation_key)
    plt.tight_layout()

    output_path = _prepare_output_path(output_path)
    if output_path is not None:
        plt.savefig(output_path, dpi=200, bbox_inches="tight")
        plt.close()
    return table


def plot_cnv_outputs(
    adata,
    cluster_keys: Sequence[str],
    output_dir: str | Path,
    spatial_key: str = "spatial",
    point_size: float = 4.0,
    heatmap_vmin: float = -0.4,
    heatmap_vmax: float = 0.4,
):
    """Save heatmap and spatial plots for one or more CNV cluster keys."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for key in cluster_keys:
        safe_key = key.replace("/", "_")
        plot_chromosome_heatmap(
            adata,
            groupby=key,
            output_path=output_dir / f"{safe_key}_heatmap.png",
            vmin=heatmap_vmin,
            vmax=heatmap_vmax,
        )
        plot_spatial(
            adata,
            color=key,
            output_path=output_dir / f"{safe_key}_spatial.png",
            spatial_key=spatial_key,
            point_size=point_size,
            title=f"Spatial {key}",
        )
