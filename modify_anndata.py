import os
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp
import pandas as pd


def data_info(adata):
    
    # Convert count matrix to dense array (if sparse)
    count_matrix = adata.X.toarray() if sp.issparse(adata.X) else adata.X
     
    # Size of the data
    print( 'Size of the data: \n',\
        len(adata.var_names), 'genes \n',\
        len(adata.obs_names), 'cells \n',\
        len(count_matrix.flatten()), 'size of count matrix \n',\
        len(count_matrix.flatten()[count_matrix.flatten() != 0]), 'positions in the count matrix are non-zero \n'\
    )

    #Check if raw data exists
    if adata.raw is None:
        print("No raw data found in the AnnData object.")
    else:
        print("Raw data found in the AnnData object.")

    gene_count_mtx = create_raw_count_mtx(adata)
    
    calculate_nonzero_average(gene_count_mtx)


def convert_to_raw(adata):
    
    # Extract the raw count matrix
    if adata.raw is None:
        raise ValueError("No raw data found in the AnnData object.")
    
    # Create a new AnnData object with raw data as count matrix
    adata_raw = ad.AnnData(X=adata.raw.X,
                                  obs=adata.obs,
                                  var=adata.var)
    
    # Copy additional annotations if needed
    adata_raw.obs = adata.obs.copy()
    adata_raw.var = adata.var.copy()
    adata_raw.uns = adata.uns.copy()
    
    return(adata_raw)



def create_raw_count_mtx(adata, output_file_path=None):
    adata = convert_to_raw(adata)
    
    count_matrix = adata.X.toarray() if sp.issparse(adata.X) else adata.X

    df = pd.DataFrame(count_matrix)
    df.columns = adata.var_names
    df.index = adata.obs_names

    if output_file_path:
        df.to_csv(output_file_path)
        print('Raw count matrix saved as', output_file_path)

    return df


def calculate_nonzero_average(gene_count_mtx):

    # Select non-zero values
    count_matrix = np.asarray(gene_count_mtx) if not isinstance(gene_count_mtx, np.ndarray) else gene_count_mtx
    non_zero_values = count_matrix.flatten()[count_matrix.flatten() != 0]
    
    # Calculate the mean of non-zero values
    mean_non_zero = np.mean(non_zero_values)
    
    print('The average read depth (excluding zeros):', mean_non_zero)
    return mean_non_zero


# Is binomial distribution the best approach to subsample the read depth?
def downsample_counts(adata, downsample_factor, output_file_path=None):
    adata_count_mtx = create_raw_count_mtx(adata)
    downsampled_matrix = adata_count_mtx.map(lambda x: np.random.binomial(n=x, p=downsample_factor))

    if output_file_path:
        downsampled_matrix.to_csv(output_file_path)
        print('Raw count matrix saved as', output_file_path)
    
    return downsampled_matrix

