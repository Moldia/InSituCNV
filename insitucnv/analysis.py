import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score, davies_bouldin_score, silhouette_score
from tqdm.auto import tqdm


def _calculate_spatial_cohesion(adata: ad.AnnData, cluster_key: str, spatial_key: str = "spatial") -> float:
    cohesion_scores = {}
    for cluster in adata.obs[cluster_key].unique():
        cluster_cells = adata.obs[cluster_key] == cluster
        if cluster_cells.sum() > 1:
            coords = adata.obsm[spatial_key][cluster_cells]
            cohesion_scores[cluster] = pdist(coords).mean()
        else:
            cohesion_scores[cluster] = np.nan

    valid_scores = [score for score in cohesion_scores.values() if not np.isnan(score)]
    if not valid_scores:
        return np.nan

    max_score = max(valid_scores) or 1.0
    normalized_scores = {key: 1 - (value / max_score) for key, value in cohesion_scores.items() if not np.isnan(value)}
    return float(np.nanmean(list(normalized_scores.values())))


def _calculate_cluster_stability(
    adata: ad.AnnData,
    resolution: float,
    obsm_key: str,
    n_subsamples: int,
    subsample_frac: float,
    n_neighbors: int = 15,
    neighbors_key: str = "cnv_neighbors",
) -> float:
    if neighbors_key not in adata.uns:
        sc.pp.neighbors(adata, use_rep=obsm_key, n_neighbors=n_neighbors, key_added=neighbors_key)

    adata_full = sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="_full",
        neighbors_key=neighbors_key,
        copy=True,
    )
    full_labels = adata_full.obs["_full"]

    ari_scores = []
    subsample_size = max(2, int(adata.n_obs * subsample_frac))
    for _ in range(n_subsamples):
        subsample_indices = np.random.choice(adata.obs_names, size=subsample_size, replace=False)
        adata_sub = adata[subsample_indices, :].copy()
        sc.pp.neighbors(adata_sub, use_rep=obsm_key, n_neighbors=n_neighbors, key_added=neighbors_key)
        sc.tl.leiden(adata_sub, resolution=resolution, key_added="_sub", neighbors_key=neighbors_key)
        ari = adjusted_rand_score(full_labels.loc[subsample_indices], adata_sub.obs["_sub"])
        ari_scores.append(ari)

    return float(np.mean(ari_scores))


def find_optimal_clustering(
    adata: ad.AnnData,
    resolutions: list[float],
    obsm_key: str = "X_cnv",
    spatial_key: str = "spatial",
    n_subsamples: int = 5,
    subsample_frac: float = 0.8,
    n_neighbors: int = 15,
    neighbors_key: str = "cnv_neighbors",
) -> pd.DataFrame:
    cnv_matrix = adata.obsm[obsm_key]
    cnv_matrix_dense = cnv_matrix.toarray() if hasattr(cnv_matrix, "toarray") else cnv_matrix

    if neighbors_key not in adata.uns:
        sc.pp.neighbors(adata, use_rep=obsm_key, n_neighbors=n_neighbors, key_added=neighbors_key)

    results = []
    for res in tqdm(resolutions, desc="Testing CNV cluster resolutions"):
        cluster_key = f"cnv_leiden_{res:g}"
        sc.tl.leiden(adata, resolution=res, key_added=cluster_key, neighbors_key=neighbors_key)
        labels = adata.obs[cluster_key]
        n_clusters = len(labels.unique())
        if n_clusters < 2:
            continue

        sil_score = silhouette_score(cnv_matrix_dense, labels)
        db_score = davies_bouldin_score(cnv_matrix_dense, labels)
        stability_score = _calculate_cluster_stability(
            adata.copy(),
            res,
            obsm_key,
            n_subsamples,
            subsample_frac,
            n_neighbors=n_neighbors,
            neighbors_key=neighbors_key,
        )
        spatial_cohesion = _calculate_spatial_cohesion(adata, cluster_key, spatial_key)

        results.append(
            {
                "resolution": res,
                "n_clusters": n_clusters,
                "silhouette_score": sil_score,
                "davies_bouldin_score": db_score,
                "stability_score": stability_score,
                "spatial_cohesion_score": spatial_cohesion,
            }
        )

    return pd.DataFrame(results)
