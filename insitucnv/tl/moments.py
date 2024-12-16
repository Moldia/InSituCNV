import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import numpy as np
from scipy.sparse import csr_matrix, issparse
from scvelo import logging as logg
from scvelo import settings
from scvelo.preprocessing.neighbors import get_connectivities, get_n_neighs, neighbors, verify_neighbors
from scvelo.preprocessing.utils import normalize_per_cell, not_yet_normalized
import infercnvpy as cnv
import matplotlib
from scipy.sparse import csr_matrix


def smooth_data_for_cnv(data, n_neighbors=20, mode="connectivities",copy=None):
    """Smooths data for CNV inference using nearest neighbor connectivities.

    Parameters
    ----------
    data : :class:`~anndata.AnnData`
        Annotated data matrix.
    n_neighbors : `int`, optional (default: 20)
        Number of neighbors to use for smoothing.
    mode : {'connectivities', 'distances'}, optional (default: 'connectivities')
        Metric to use for smoothing computations.

    Returns
    -------
    None
        Modifies the input AnnData object in place by adding smoothed data to `adata.layers['M']`.
    """

    # Ensure neighbor graph is computed if required
    if n_neighbors > get_n_neighs(data):
        verify_neighbors(data)

    # Compute smoothing based on the specified mode
    connectivities = get_connectivities(data, mode, n_neighbors=n_neighbors, recurse_neighbors=False)
    data.layers["M"] = (
        csr_matrix.dot(connectivities, csr_matrix(data.layers["raw"])).astype(np.float32).toarray()
    )
    
    return adata if copy else None
