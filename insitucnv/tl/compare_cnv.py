import pandas as pd
import numpy as np
from scipy.sparse import issparse
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt



def compute_avg_cnv_profiles(adata, group_by, label):
    '''
    Function to compute the average CNV values of a dataset run by inferCNV. 
    
    Parameters:
    - adata (AnnData object): containing adata.obsm['X_cnv'] and adata.uns['cnv']['chr_pos'] after running inferCNV.
    - group_by (str): column in adata.obs to group by when calculating average CNVs. 
    - label (str): prefix for cluster labels, e.g. 'CosMx' --> CosMx_0, CosMx_1, etc. 
    
    Returns:
    - dict of DataFrames: each containing the average CNV values for each genomic bin per cluster
    '''
    
    # Ensure input is valid
    if 'X_cnv' not in adata.obsm:
        raise ValueError("Missing 'X_cnv' in adata.obsm")
    if 'cnv' not in adata.uns or 'chr_pos' not in adata.uns['cnv']:
        raise ValueError("Missing 'cnv' or 'chr_pos' in adata.uns")

    # Get unique cluster values
    adata_clusters = adata.obs[group_by].unique()

    cluster_CNV_profiles = {}

    for cluster in adata_clusters:
        # Subset the AnnData object for the cluster
        adata_cluster = adata[adata.obs[group_by] == cluster]

        X_cnv = adata_cluster.obsm['X_cnv']
        chr_pos = adata.uns['cnv']['chr_pos']  # use global chr_pos, not per cluster

        num_bins = X_cnv.shape[1]
        bin_to_chr = pd.Series(index=np.arange(num_bins), dtype="object")

        # Sort chromosomes by their start index
        sorted_chr = sorted(chr_pos.items(), key=lambda x: x[1])
        for i, (chrom, start) in enumerate(sorted_chr):
            end = sorted_chr[i + 1][1] if i + 1 < len(sorted_chr) else num_bins
            bin_to_chr[start:end] = chrom

        # Convert to dense array if sparse
        X_dense = X_cnv.toarray() if issparse(X_cnv) else X_cnv

        # Create DataFrame for cluster
        df = pd.DataFrame({
            'bin_index': np.arange(num_bins),
            'chromosome': bin_to_chr.values,
            'mean_cnv': X_dense.mean(axis=0),
            'std_cnv': X_dense.std(axis=0),
            'source': f'{label}_{cluster}'
        })

        # Store in dictionary
        cluster_CNV_profiles[f'{label}_{cluster}'] = df

    return cluster_CNV_profiles



def compute_cnv_similarity_matrix(profiles_a, profiles_b, sort_bins=True):
    '''
    Compute a cosine similarity matrix between two sets of CNV cluster profiles.

    Parameters:
    - profiles_a (dict): Output from summarize_cnv for dataset A (e.g., snPATHOseq)
    - profiles_b (dict): Output from summarize_cnv for dataset B (e.g., CosMx)
    - sort_bins (bool): If True, sorts DataFrames by 'bin_index' to ensure alignment

    Returns:
    - similarity_matrix (pd.DataFrame): Cosine similarity matrix between clusters
    '''

    # Sort each DataFrame by bin_index to ensure vectors are aligned
    if sort_bins:
        for df in profiles_a.values():
            df.sort_values(by='bin_index', inplace=True)
        for df in profiles_b.values():
            df.sort_values(by='bin_index', inplace=True)

    # Initialize empty similarity matrix
    similarity_matrix = pd.DataFrame(
        index=profiles_a.keys(),
        columns=profiles_b.keys(),
        dtype=float
    )

    # Compute cosine similarities
    for label_a, df_a in profiles_a.items():
        vec_a = df_a['mean_cnv'].values.reshape(1, -1)

        for label_b, df_b in profiles_b.items():
            vec_b = df_b['mean_cnv'].values.reshape(1, -1)

            similarity = cosine_similarity(vec_a, vec_b)[0, 0]
            similarity_matrix.loc[label_a, label_b] = similarity

    return similarity_matrix





# Helper function to compute averaged CNV from multiple profiles
def average_profiles(cluster_list, profile_dict):
    '''
    Computes the average CNV profile across a list of clusters.

    Parameters:
    - cluster_list (list of str): List of cluster labels to include in the averaging. 
      Each label must correspond to a key in `profile_dict`.
    - profile_dict (dict of pd.DataFrame): Dictionary of CNV profiles per cluster,
      typically the output from `compute_avg_cnv_profiles`. Each DataFrame must
      contain 'bin_index', 'chromosome', and 'mean_cnv' columns.

    Returns:
    - avg_df (pd.DataFrame): Containing the average CNV profile across the selected clusters. Includes 'bin_index', 'chromosome', and 
      the averaged 'mean_cnv' values.
    '''
    
    dfs = [profile_dict[label].sort_values('bin_index').reset_index(drop=True) for label in cluster_list]
    # Sanity check: ensure same bin structure
    assert all(np.array_equal(dfs[0]['bin_index'].values, df['bin_index'].values) for df in dfs)

    # Stack mean CNV values and compute average
    mean_cnv_matrix = np.stack([df['mean_cnv'].values for df in dfs], axis=0)
    avg_cnv = mean_cnv_matrix.mean(axis=0)

    # Use bin_index and chromosome from the first df
    avg_df = dfs[0][['bin_index', 'chromosome']].copy()
    avg_df['mean_cnv'] = avg_cnv

    return avg_df


    

def plot_avg_cnv_comparison(df1, df2, chr_pos, label1='Dataset A', label2='Dataset B', title='Average CNV Comparison', save_path=None):
    '''
    Plot a side-by-side comparison of two averaged CNV profiles, highlighting gains, losses, and regions where the CNV direction differs between the two datasets.

    Parameters:
    - df1 (pd.DataFrame): CNV profile for the first dataset. Must contain columns 'bin_index', 'mean_cnv', and 'chromosome'.
    - df2 (pd.DataFrame): CNV profile for the second dataset. Must contain columns 'bin_index', 'mean_cnv', and 'chromosome'.
    - chr_pos (dict): Dictionary mapping chromosome names to their start bin index, typically from `adata.uns['cnv']['chr_pos']`.
    - label1 (str): Label for the first dataset (used in the plot legend).
    - label2 (str): Label for the second dataset (used in the plot legend).
    - title (str): Title of the plot.
    - save_path (str or None): If provided, the plot will be saved to this file path (e.g., 'output/cnv_plot.pdf'). If None, the plot is shown interactively.

    Returns:
    - None. Displays or saves the plot.
    '''
    sorted_chr_items = sorted(chr_pos.items(), key=lambda x: x[1])

    x_vals = df1['bin_index'].values

    y1 = df1['mean_cnv'].values
    y2 = df2['mean_cnv'].values

    plt.figure(figsize=(20, 6))
    plt.plot(x_vals, y1, color='gray', linestyle=':', label=label1, linewidth=2)
    plt.plot(x_vals, y2, color='gray', linestyle='-', label=label2, linewidth=2)

    # Highlight gains and losses
    plt.fill_between(x_vals, 0, y1, where=(y1 > 0), color='red', alpha=0.3)
    plt.fill_between(x_vals, 0, y1, where=(y1 < 0), color='blue', alpha=0.3)
    plt.fill_between(x_vals, 0, y2, where=(y2 > 0), color='red', alpha=0.3)
    plt.fill_between(x_vals, 0, y2, where=(y2 < 0), color='blue', alpha=0.3)

    # Highlight discordant CNV directions
    opposite_sign = (np.sign(y1) != np.sign(y2)) & (np.sign(y1) != 0) & (np.sign(y2) != 0)
    plt.fill_between(x_vals, np.minimum(y1, y2), np.maximum(y1, y2),
                     where=opposite_sign, color='yellow', alpha=1)

    # Draw chromosome boundaries
    for chrom, start in sorted_chr_items:
        plt.axvline(start, color='gray', linestyle='-', linewidth=0.5)

    # Add chromosome labels
    chrom_labels = []
    label_positions = []
    for j, (chrom, start) in enumerate(sorted_chr_items):
        end = sorted_chr_items[j + 1][1] if j + 1 < len(sorted_chr_items) else x_vals.max()
        center = (start + end) // 2
        chrom_labels.append(chrom)
        label_positions.append(center)

    plt.xticks(label_positions, chrom_labels, rotation=90)
    plt.xlabel('Genomic bin')
    plt.ylabel('Mean CNV')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.xlim(x_vals.min(), x_vals.max())


    if save_path:
        plt.savefig(save_path)
    plt.show()