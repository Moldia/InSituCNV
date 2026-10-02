from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix, identity


def smooth_data_for_cnv(
    data,
    layer_w_norm_counts: str | None = None,
    n_neighbors: int = 20,
    mode: str = "connectivities",
    output_layer: str = "M",
    neighbors_key: str = "neighbors",
    copy: bool = False,
):
    """Smooth expression values over a nearest-neighbor graph before CNV inference.

    This is the reusable version of the notebook step:
    normalized counts -> graph smoothing -> smoothed counts in ``adata.layers["M"]``.

    Each cell is replaced by the unweighted mean of itself and its neighbors in the
    graph. The graph must have been built with exactly ``n_neighbors`` (e.g.
    ``sc.pp.neighbors(adata, n_neighbors=n_neighbors)``); a mismatch raises an error
    instead of silently smoothing over a different neighborhood.

    Parameters
    ----------
    data
        AnnData object with a precomputed neighbor graph.
    layer_w_norm_counts
        Layer containing normalized, non-log-transformed counts. If ``None``,
        ``adata.X`` is smoothed.
    n_neighbors
        Number of nearest neighbors the graph was built with. Follows the scanpy
        convention: the count includes the cell itself. With
        ``mode="connectivities"`` the graph is symmetric, so cells that are
        neighbors of many others average over more than ``n_neighbors`` cells.
    mode
        Graph to use, ``"connectivities"`` (symmetric) or ``"distances"`` (the
        plain kNN graph: each cell plus its ``n_neighbors - 1`` nearest).
    output_layer
        Layer where smoothed counts are stored. The manuscript notebooks used
        ``"M"``, which remains the default for compatibility.
    neighbors_key
        ``key_added`` of the scanpy neighbor graph to use (``adata.uns[neighbors_key]``).
    copy
        Return a modified copy instead of editing ``data``.

    Returns
    -------
    AnnData
        The modified AnnData object.
    """
    if mode not in ("connectivities", "distances"):
        raise ValueError(f"mode must be 'connectivities' or 'distances', got {mode!r}.")

    adata = data.copy() if copy else data

    graph_info = adata.uns.get(neighbors_key)
    if graph_info is None:
        raise KeyError(
            f"No neighbor graph found in adata.uns['{neighbors_key}']. Run "
            f"sc.pp.neighbors(adata, n_neighbors={n_neighbors}) first."
        )
    graph_n_neighbors = graph_info.get("params", {}).get("n_neighbors")
    if graph_n_neighbors != n_neighbors:
        raise ValueError(
            f"The neighbor graph in adata.uns['{neighbors_key}'] was built with "
            f"n_neighbors={graph_n_neighbors}, but smoothing requested n_neighbors={n_neighbors}. "
            f"Rebuild it with sc.pp.neighbors(adata, n_neighbors={n_neighbors}), or use "
            "insitucnv.tl.prepare_cnv_input, which builds a matching graph."
        )
    graph = adata.obsp[graph_info.get(f"{mode}_key", mode)]

    # Uniform average over the cell and its graph neighbors (same as scvelo's
    # get_connectivities without neighbor trimming).
    weights = csr_matrix(((graph > 0) + identity(adata.n_obs, format="csr")) > 0, dtype=np.float32)
    weights = csr_matrix(weights.multiply(1.0 / weights.sum(axis=1)))

    if layer_w_norm_counts is None:
        matrix = adata.X
    else:
        if layer_w_norm_counts not in adata.layers:
            raise KeyError(
                f"Layer '{layer_w_norm_counts}' not found in adata.layers. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        matrix = adata.layers[layer_w_norm_counts]

    adata.layers[output_layer] = csr_matrix.dot(weights, csr_matrix(matrix)).astype(np.float32).toarray()
    return adata
