This folder contains the code used to generate the outputs presented in Figure 1 and some of the supplementary figures of the manuscript "*In Situ Inference of Copy Number Variations in Image-Based Spatial Transcriptomics*".

<img src="https://github.com/user-attachments/assets/356d4fb4-7884-4deb-8128-b4a100a8b6a2" width="700">

**a**, CNV heatmap of the colorectal carcinoma CosMx sample, showing gains (red) and losses (blue), including
epithelial cells grouped by CNV cluster (0-3), and non-epithelial cells (N). Malignant subclones are marked in the
dendrogram. **b**, Spatial map of epithelial CNV clusters presented in a. **c**, Histopathological annotations of the CosMx
section, identifying healthy reference-like mucosa (REF), tubulovillous adenomas (TVA) stage 1-3, and cancerous
colorectal carcinoma (CRC) regions. **d**, CNV cluster composition across histopathological compartments. **e**,
Genome-wide comparison of average CNV profiles for the malignant CosMx (solid line) and snPATHO-seq (dotted
line) cells, highlighting common gains (red), losses (blue) and disagreement between the datasets (yellow). **f**, CNV
heatmap of the HGOC sample, illustrating inferred gains (red) and losses (blue). **g**, Map of cell classes in the HGOC
dataset. **h**, Map of epithelial CNV clusters defined in f.. **i**, Relative frequencies of T cells, B cells, and fibroblasts
plotted as a function of their distance to the nearest malignant subclone, as defined in h.

## Folder structure

### CosMx_CRC 
Includes the analysis performed on the CosMx and snPATHO-seq datasets from colorectal tumor sample 221, avaliable in Crowell et al., "*Tracing colorectal malignancy transformation from cell to tissue scale*", bioRxiv (2025) https://doi.org/10.1101/2025.06.23.660674.
- 1_R_convert_from_SCE_extract_metadata
- 2_Prepare_data
- 3_inferCNV
- 4_Compare

### Xenium_prime_ovarian 
Contains analysis of the publicly available high-grade ovarian cancer Xenium Prime dataset from 10x Genomics: https://www.10xgenomics.com/datasets/xenium-prime-fresh-frozen-human-ovary.

- 1_Cell_typing_Ovarian.ipynb
- 2_infercnv_Ovarian_moments_correct.ipynb
- 3_infercnv_ovarian_moments_ACROSS_winsize_n_read_depth.ipynb
- 4_spatial_stats_Ovarian_moments.ipynb

### Xenium_prime_lymph_node
Contains analysis of the publicly available lymph node Xenium Prime dataset from 10x Genomics: https://www.10xgenomics.com/datasets/preview-data-xenium-prime-gene-expression.
- 1_Cell_typing_lymphnode.ipynb
- 2_infercnv_lymphnode_moments.ipynb
