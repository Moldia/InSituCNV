import numpy as np
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import warnings

# Suppress warnings to avoid cluttered output
warnings.simplefilter("ignore")

# Set default figure parameters for plotting
sc.settings.set_figure_params(figsize=(5, 5))



def inferCNV_for_simulated_data(path, data_name, window_size):
    """
    Perform CNV (Copy Number Variation) inference on a simulated dataset using infercnvpy.

    Parameters:
    - path (str): The directory path where the dataset is stored.
    - data_name (str): The name of the dataset (without file extension).
    - window_size (int): The size of the genomic window for CNV inference.
    - step (int): The step size for sliding window CNV calculation.

    Outputs:
    - Saves CNV heatmap and UMAP visualizations.
    - Saves the processed AnnData object with inferred CNVs.
    """

    # Load the dataset from the specified path
    adata = sc.read_h5ad(path + '/' + data_name + '.h5ad')
    
    # Extract the simulated CNV layer and assign it to the main matrix (X)
    layer = 'CNV_simulated'
    adata.X = adata.layers[layer].copy()

    # Perform CNV inference using infercnvpy
    cnv.tl.infercnv(
        adata,
        reference_key="cell_type",  # Column in adata.obs containing reference cell types
        reference_cat=['ciliated cell', 'basal cell'],  # List of reference cell types
        window_size=window_size,
        step=1,
        calculate_gene_values=True,  # Compute CNV values for genes
    )

    # Perform dimensionality reduction and clustering
    cnv.tl.pca(adata)  # Compute PCA for feature reduction
    cnv.pp.neighbors(adata)  # Compute nearest neighbors graph
    cnv.tl.leiden(adata)  # Perform clustering using the Leiden algorithm

    # Adjust clustering resolution until exactly four clusters are obtained
    resol = 0.7
    cnv.tl.leiden(adata, resolution=resol, key_added='cnv_leiden')

    while len(np.unique(adata.obs['cnv_leiden'])) != 4:
        print('Adjusting the resolution of the clustering...')
        
        if len(np.unique(adata.obs['cnv_leiden'])) < 4:
            resol += 0.005  # Increase resolution if there are fewer than 4 clusters
            print(f'Increasing resolution to {resol}')
        elif len(np.unique(adata.obs['cnv_leiden'])) > 4:
            resol -= 0.005  # Decrease resolution if there are more than 4 clusters
            print(f'Reducing resolution to {resol}')

        cnv.tl.leiden(adata, resolution=resol, key_added='cnv_leiden')

    # Store final clustering results with the adjusted resolution
    cnv.tl.leiden(adata, resolution=resol, key_added=f'cnv_leiden_res{round(resol, 2)}')

    # Generate and save a CNV heatmap
    cnv.pl.chromosome_heatmap(
        adata,
        groupby="cnv_leiden",  # Group cells by inferred CNV clusters
        dendrogram=False,  # Do not plot a dendrogram
        vmin=-0.4, vmax=0.4,  # Set color range for CNV values
        save=f'_{data_name}.png'
    )
    print('CNV heatmap saved!')

    # Compute and save UMAP visualizations
    cnv.tl.umap(adata)  # Compute UMAP embedding
    cnv.tl.cnv_score(adata, groupby='cnv_leiden')  # Compute CNV scores per cluster

    cnv.pl.umap(
        adata,
        color=['simulated_subclone', 'cnv_leiden', 'cnv_score'],  # Colors for UMAP plots
        save=f'_{data_name}.png'
    )
    print('CNV UMAP saved!')

    # Save the processed AnnData object with inferred CNVs
    new_path = path + '/' + data_name + '_CNVinf.h5ad'
    adata.write(new_path, compression='gzip')

    print(f'CNV inference for the {data_name} dataset has been completed. The results have been saved in: {new_path}')




def save_figures(path, data_name):
        '''
    Loads a CNV-inferred AnnData object and generates visualizations of CNV heatmaps and UMAPs,
    saving them as PDF files.

    Parameters:
        path (str): Directory path where the AnnData object is stored.
        data_name (str): Base name of the AnnData file (excluding the '_CNVinf.h5ad' suffix).

    Returns:
        None
        '''

    # Load the dataset from the specified path
    adata = sc.read_h5ad(path + '/' + data_name + '_CNVinf.h5ad')

    cnv_leiden_colors = ["#df9a57","#fc7a57","#fcd757","#a22020"]
    simulated_subclone_colors = ["#bce784","#5dd39e","#348aa7","#525174"]
    
    # Assign custom colors to the respective categories
    adata.uns["cnv_leiden_colors"] = cnv_leiden_colors
    adata.uns["simulated_subclone_colors"] = simulated_subclone_colors
    
    # Generate and save a CNV heatmap
    cnv.pl.chromosome_heatmap(
        adata,
        groupby="cnv_leiden",  # Group cells by inferred CNV clusters
        dendrogram=False,  # Do not plot a dendrogram
        vmin=-0.4, vmax=0.4,  # Set color range for CNV values
        save=f'_{data_name}.pdf'
    )
    print('CNV heatmap saved!')

    cnv.pl.umap(
        adata,
        color=[ 'cnv_leiden','simulated_subclone'],  # Colors for UMAP plots
        save=f'_{data_name}.pdf'
    )
    print('CNV UMAP saved!')
