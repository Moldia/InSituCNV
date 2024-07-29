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
        len(count_matrix.flatten()[count_matrix.flatten() != 0]), 'positive positions in the count matrix \n'\
    )

    #Check if raw data exists
    if adata.raw is None:
        print("No raw data found in the AnnData object. \n")
    else:
        print("Raw data found in the AnnData object. \n")
        adata.layers['raw'] = adata.raw.X.toarray()

    gene_count_mtx = create_raw_count_mtx(adata)

    print('Average read depths:')
    calculate_read_depths(adata)


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



def create_raw_count_mtx(adata, layer_name = None, output_file_path=None):
    if layer_name:
        if layer_name not in adata.layers:
            raise ValueError(f"No layer named '{layer_name}' found in the AnnData object.")

        count_matrix = adata.layers[layer_name]
    
        if sp.issparse(count_matrix):
            count_matrix = count_matrix.toarray()
        
    else: 
        adata_raw = convert_to_raw(adata)
        count_matrix = adata_raw.X.toarray() if sp.issparse(adata_raw.X) else adata_raw.X
    
    df = pd.DataFrame(count_matrix)
    df.columns = adata.var_names
    df.index = adata.obs_names

    if output_file_path:
        if not os.path.exists(output_file_path):
            df.to_csv(output_file_path)
            print('Raw count matrix saved as', output_file_path)
        else:
            print('Raw count matrix already exists!')
        
    return df





# Is binomial distribution the best approach to subsample the read depth?
def downsample_counts(adata, downsample_factor, output_file_path=None):

    if output_file_path:
        if not os.path.exists(output_file_path):
            adata_count_mtx = create_raw_count_mtx(adata)
            downsampled_matrix = adata_count_mtx.map(lambda x: np.random.binomial(n=x, p=downsample_factor))
            
            downsampled_matrix.to_csv(output_file_path)
            print('Raw count matrix saved as', output_file_path)

            return downsampled_matrix
        else:
            print('Raw count matrix already exists!')
    
    else:
        adata_count_mtx = create_raw_count_mtx(adata)
        downsampled_matrix = adata_count_mtx.map(lambda x: np.random.binomial(n=x, p=downsample_factor))
        
        return downsampled_matrix




def calculate_read_depths(adata):
    
    # Obtain the indices of the original positive values
    positive_indices = np.nonzero(adata.layers['raw'])

    # Calculate the average value at these indices across all layers
    average_values = {}
    for layer_name in adata.layers.keys():
        layer_matrix = adata.layers[layer_name]
        positive_values = layer_matrix[positive_indices]
        average_values[layer_name] = positive_values.mean()
        
        # Print the average values for each layer
        print(layer_name,':', round(positive_values.mean(), 2))
    
    return average_values






def test_exists(path):
    if not os.path.exists(path):
        print('the path did not exist!')
    else:
        print('the path did exist!')
