from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix
from scvelo.preprocessing.neighbors import get_connectivities, get_n_neighs, verify_neighbors


def smooth_data_for_cnv(
    data,
    layer_w_norm_counts: str | None = None,
    n_neighbors: int = 20,
    mode: str = "connectivities",
    output_layer: str = "M",
    copy: bool = False,
):
    """Smooth expression values over a nearest-neighbor graph before CNV inference.

    This is the reusable version of the notebook step:
    normalized counts -> graph smoothing -> smoothed counts in ``adata.layers["M"]``.

    Parameters
    ----------
    data
        AnnData object with a precomputed neighbor graph.
    layer_w_norm_counts
        Layer containing normalized, non-log-transformed counts. If ``None``,
        ``adata.X`` is smoothed.
    n_neighbors
        Number of nearest neighbors to use from the graph.
    mode
        Graph weights to use, usually ``"connectivities"`` or ``"distances"``.
    output_layer
        Layer where smoothed counts are stored. The manuscript notebooks used
        ``"M"``, which remains the default for compatibility.
    copy
        Return a modified copy instead of editing ``data``.

    Returns
    -------
    AnnData
        The modified AnnData object.
    """
    adata = data.copy() if copy else data

    if n_neighbors > get_n_neighs(adata):
        verify_neighbors(adata)

    connectivities = get_connectivities(adata, mode, n_neighbors=n_neighbors, recurse_neighbors=False)
    if layer_w_norm_counts is None:
        matrix = adata.X
    else:
        if layer_w_norm_counts not in adata.layers:
            raise KeyError(
                f"Layer '{layer_w_norm_counts}' not found in adata.layers. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        matrix = adata.layers[layer_w_norm_counts]

    adata.layers[output_layer] = csr_matrix.dot(connectivities, csr_matrix(matrix)).astype(np.float32).toarray()
    return adata
