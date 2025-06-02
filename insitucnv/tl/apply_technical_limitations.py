import scanpy as sc
import numpy as np

def subsample_counts(adata, fraction, destination, data_name):
    """
    Subsamples the total counts in an AnnData object to simulate lower sequencing depth,
    then normalizes, log-transforms, and saves the result to disk.

    Parameters:
        adata (AnnData object): AnnData object containing gene expression data to be subsampled.
        fraction (nbr): A value between 0 and 1 indicating the proportion of total counts to retain.
                        For example, 0.5 will retain 50% of the total counts.
        destination (str): Directory path where the resulting AnnData object will be saved.
        data_name (str): Filename (without extension) for the saved AnnData object.

    Returns:
        None

    """
    # Creates a copy of the input `adata` to preserve the original data
    new_adata = adata.copy()
    layer_name = 'CNV_simulated'

    # Validates that `fraction` is within the [0, 1] range
    if fraction > 1 or fraction < 0:
        raise ValueError(f"`fraction` needs to be within [0, 1], not {fraction}")

    if fraction != 1:
        total_counts = int(fraction * adata.X.sum()) 
        downsampled_adata = sc.pp.downsample_counts(
            adata, 
            total_counts=total_counts, 
            random_state=42, 
            copy=True
        )

        new_adata.layers[layer_name] = downsampled_adata.X.copy()

    # If `fraction` == 1, the original counts are used unchanged
    else: 
        new_adata.layers[layer_name] = new_adata.X.copy()

    # Normalize and log-transform the downsampled data
    new_adata.layers[f"{layer_name}_raw"] = new_adata.layers[layer_name].copy()
    sc.pp.normalize_total(new_adata, layer=layer_name)
    sc.pp.log1p(new_adata, layer=layer_name)

    new_adata.X = new_adata.layers[layer_name].copy()

    # Save the AnnData object to the specified location
    new_adata.write(destination + '/' + data_name + '.h5ad', compression='gzip')
    print(f'subsampled data ({fraction*100:.0f}% counts) saved as {destination}/{data_name}.h5ad')



def subsample_genes(adata, gene_panel_size, destination, data_name):
    """
    Randomly subsamples genes from an AnnData object to simulate reduced gene panels,
    and saves the resulting object to disk.

    Parameters:
        adata (AnnData object): AnnData object containing gene expression data to be subsampled.
        gene_panel_size (nbr or 'all'): The number of genes to retain. If 'all', no subsampling is performed.
        destination (str): Directory path where the resulting AnnData object will be saved.
        data_name (str): Filename (without extension) for the saved AnnData object.

    Returns:
        None
    """
    
    # Set random seed for reproducibility
    random_state = 42
    np.random.seed(random_state)

    # Calculate the number of genes to retain
    old_n_vars = adata.n_vars
    
    if gene_panel_size == 'all':
        gene_panel_size = old_n_vars
        
    new_n_vars =  gene_panel_size
    var_indices = np.random.choice(old_n_vars, size=new_n_vars, replace=False)

    # Subsample the data to include only selected genes
    subsampled_adata = adata[:, var_indices].copy()
                
    # Save the AnnData object to the specified location
    subsampled_adata.write(destination + '/' + data_name + '.h5ad', compression='gzip')
    print(f'subsampled data (with panel size:{new_n_vars} genes) saved as {destination}/{data_name}.h5ad')
