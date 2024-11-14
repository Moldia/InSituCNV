import scanpy as sc
import pandas as pd
import infercnvpy as cnv # Assuming the cnv package has the necessary infercnv functionality

# Step 1: Run inferCNV
def step1_infercnv(adata, layer, window_size, reference_key, reference_cat):
    """
    Perform CNV inference on the given AnnData object.
    
    Parameters:
        adata: AnnData object
        layer: Data layer to use ('raw', 'log_counts', etc.)
        window_size: Window size for CNV analysis
        reference_key: Key used for reference (e.g., 'leiden')
        reference_cat: Categories within the reference (e.g., [3, 5, 17])
    
    Returns:
        Updated AnnData object with inferred CNV information.
    """
    print(f"Running inferCNV with layer={layer}, window_size={window_size}, reference_key={reference_key}, reference_cat={reference_cat}")
    
    cnv.tl.infercnv(
        adata,
        layer=layer,
        window_size=window_size,
        reference_key=reference_key,
        reference_cat=reference_cat
    )
    
    return adata


# Step 2: CNV Leiden Clustering
def step2_cnv_leiden_clustering(adata, cnv_leiden_resolutions):
    """
    Perform Leiden clustering based on CNV at multiple resolutions.
    
    Parameters:
        adata: AnnData object
        cnv_leiden_resolutions: List of resolutions for Leiden clustering
    
    Returns:
        Updated AnnData object with CNV Leiden clusters.
    """
    print(f"Running CNV Leiden clustering with resolutions: {cnv_leiden_resolutions}")
    heatmap_groupby_additions = []
    
    for res in cnv_leiden_resolutions:
        key_added = f'cnv_leiden_res{res}'
        heatmap_groupby_additions.append(key_added)
        
        print(f"Running Leiden clustering with resolution={res}")
        cnv.tl.leiden(adata, resolution=res, flavor="igraph", n_iterations=2, key_added=key_added)
    
    return adata, heatmap_groupby_additions


# Step 3: Plot and save CNV heatmap
def step3_heatmap(adata, heatmap_groupby, save_heatmap):
    """
    Plot and save heatmap for CNV data.
    
    Parameters:
        adata: AnnData object
        heatmap_groupby: Grouping to use for the heatmap ('leiden', 'cell_type', etc.)
        save_heatmap: Filename to save the heatmap
    
    Returns:
        None
    """
    print(f"Generating heatmap for groupby={heatmap_groupby}, saving as {save_heatmap}")
    sc.tl.dendrogram(adata, groupby=heatmap_groupby)
    cnv.pl.chromosome_heatmap(adata, groupby=heatmap_groupby, dendrogram=True, save=save_heatmap)



## Nest steps: add the spatial image and plot it according to the leiden things generated!!

def step4_spatialplot(adata, spot_size, heatmap_groupby, save_spatialplot):
    """
    Plot and save spatial image.
    
    Parameters:
        adata: AnnData object
        spot_size: Spot size
        heatmap_groupby: Use the same grouping as for the heatmap ('leiden', 'cell_type', etc.)
        save_spatialplot: Filename to save the heatmap
    
    Returns:
        None
    """

    print(f"Generating spatial plot for groupby={heatmap_groupby}, saving as {save_spatialplot}")
    sc.settings.set_figure_params(figsize=(20, 20))
    sc.pl.spatial(adata, spot_size = spot_size, color=heatmap_groupby, save=save_spatialplot)




# Main pipeline function
def cnv_pipeline(
    adata, sample_name, layer, window_size, reference_key, reference_cat,
    cnv_leiden_resolutions, heatmap_groupby_list, spot_size
):
    """
    Full CNV analysis pipeline: inferCNV, Leiden clustering, and heatmap generation.
    
    Parameters:
        adata: AnnData object
        sample_name: Name of the sample (for file naming)
        layer: Data layer to use ('raw', 'log_counts', etc.)
        window_size: Window size for CNV analysis
        reference_key: Key used for reference (e.g., 'leiden')
        reference_cat: Categories within the reference (e.g., [3, 5, 17])
        cnv_leiden_resolutions: List of resolutions for Leiden clustering
        heatmap_groupby_list: List of groupings for heatmap ('leiden', 'cell_type', etc.)
        spot_size: Spot size
    
    Returns:
        None
    """
    # Step 1: Run inferCNV
    adata = step1_infercnv(adata, layer, window_size, reference_key, reference_cat)

    # Step 2: Run PCA and neighbors, then CNV Leiden clustering
    print("Running PCA and computing neighbors...")
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    
    adata, heatmap_groupby_additions = step2_cnv_leiden_clustering(adata, cnv_leiden_resolutions)
    
    # Automatically include CNV Leiden clusterings in heatmap groupby
    heatmap_groupby_list += heatmap_groupby_additions

    # Step 3: Generate and save heatmaps for each groupby condition
    for groupby in heatmap_groupby_list:
        save_heatmap = f"_sample_{sample_name}_{layer}_w{window_size}_ref-{reference_key}-{'-'.join(map(str, reference_cat))}_groupby-{groupby}.png"
        save_spatialplot = f"_sample_{sample_name}_{layer}_w{window_size}_ref-{reference_key}-{'-'.join(map(str, reference_cat))}_groupby-{groupby}.png"
        step3_heatmap(adata, groupby, save_heatmap)
        step4_spatialplot(adata, spot_size, groupby, save_spatialplot)

    print("Pipeline complete.")



# Example run of the pipeline
'''
adata = sc.read_h5ad("brain_tumor_data.h5ad")  # Load your AnnData object here

# Parameters
sample_name = "brainTumor10K"
layer = "log_counts"  # Or "raw"
window_size = 150
reference_key = "cell_type"
reference_cat = ["Astrocytes", "Oligodendrocytes"]
cnv_leiden_resolutions = [0.3, 0.6, 1.0]
heatmap_groupby_list = ["cell_type", "leiden"]
spot_size = 20

# Run the pipeline
cnv_pipeline(adata, sample_name, layer, window_size, reference_key, reference_cat, cnv_leiden_resolutions, heatmap_groupby_list, spot_size)
'''
