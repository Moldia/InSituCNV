import os
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp
import pandas as pd


def data_info(adata):
    '''This function writes out an overview of the data information, such as the number of genes 
    and cells present in the data, the size of the count matrix (cells x genes), and the number of 
    positive positions, i.e. where there have been observed a count. 

    Params:
        adata - an AnnData object
    '''

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



def convert_to_raw(adata):
    '''Converts an AnnData object to the raw version (if such exists), where the adata.X will contain the raw counts, 
    instead of the processed.

    Params:
        adata: an AnnData object
    Output:
        adata_raw: an AnnData object where the adata.X contains the raw counts

    '''
    
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
    '''
    
    '''
    
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







## Calculate read depths

def calculate_read_depths(adata):
    
    # Obtain the indices of the original positive values
    positive_indices = np.nonzero(adata.layers['raw'])

    # Calculate the average value at these indices across all layers
    average_values = pd.DataFrame()
    for layer_name in adata.layers.keys():
        layer_matrix = adata.layers[layer_name]
        positive_values = layer_matrix[positive_indices]
        average_values[layer_name] = positive_values.mean()
    
    return average_values

def calculate_read_depths_v2(adata, percentiles):
    # Obtain the indices of the original positive values
    positive_indices = np.nonzero(adata.layers['raw'])
    
    # Initialize a dictionary to store results
    results = {}
    
    # Loop through each layer in adata
    for layer in adata.layers.keys():
        layer_matrix = adata.layers[layer]
        
        # Average calculation
        positive_values = layer_matrix[positive_indices]
        average = positive_values.mean()
        
        # Initialize dictionary to store percentiles for this layer
        layer_percentiles = {'average': round(average, 2)}
        
        # Calculate percentiles
        for percentile in percentiles:
            gene_percentile_list = []
            for gene_index in range(layer_matrix.shape[1]):
                # Extract the counts for the current gene
                if isinstance(layer_matrix, np.ndarray):
                    gene_counts = layer_matrix[:, gene_index]
                else:
                    gene_counts = layer_matrix[:, gene_index].toarray().flatten()  # Ensure it's a flat array if sparse
                
                # Calculate the specified percentile for the current gene
                gene_percentile_value = np.percentile(gene_counts, percentile)
                gene_percentile_list.append(gene_percentile_value)
                
            # Store the average percentile value across all genes
            average_percentile = np.mean(gene_percentile_list)
            layer_percentiles[f'percentile_{percentile}'] = round(average_percentile, 2)
        
        # Store the results for this layer
        results[layer] = layer_percentiles
    
    # Convert the results dictionary to a DataFrame
    results_df = pd.DataFrame(results).T
    
    return results_df









def add_genomic_positions(adata):
    ''' Adds gene positions to the AnnData object. The adata.var_names have to be EnsmblIDs. The data is taken from the Ensmbl BioMart human dataset GRCh38.p14. 

    param:
        adata (AnnData object): the annotated dataset where the gene postions should be added
    output: 
        adata (AnnData object): the annotated dataset where the genes are added
 
    '''

    
    # Load the gene positions data
    gene_file = "/home/augusta/SSS_mount/insituCNV/InSituCNV/Ensmbl_BioMart_gene_info.txt"
    gene_positions_df = pd.read_csv(gene_file)

    # Create a dictionary for quick lookup
    gene_dict = gene_positions_df.set_index("Gene stable ID")[["Chromosome/scaffold name", "Gene start (bp)", "Gene end (bp)"]].to_dict("index")

    # Format chromosome names
    format_chromosome = lambda x: f"chr{x}"

    # Initialize counters
    genes_identified = 0

    # Iterate over the genes in adata
    for gene_id in adata.var_names:
        if gene_id in gene_dict:
            info = gene_dict[gene_id]
            adata.var.loc[gene_id, ["chromosome", "start", "end"]] = [format_chromosome(info["Chromosome/scaffold name"]), info["Gene start (bp)"], info["Gene end (bp)"]]
            genes_identified += 1

    # Print summary of identified genes
    print(f"{genes_identified} gene positions identified \n{adata.shape[1] - genes_identified} were not found")
    return adata
















