The following code was used to generate the outputs presented in Figure 1 and some of the supplementary figures of the manuscript "*In Situ Inference of Copy Number Variations in Image-Based Spatial Transcriptomics*".

![Figure 1 (1)-1](https://github.com/user-attachments/assets/356d4fb4-7884-4deb-8128-b4a100a8b6a2)
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

Figure 1 - CosMx and snPATHO-seq (Colorectal cancer)
- 1_R_convert_from_SCE_extract_metadata
- 2_Prepare_data
- 3_inferCNV
- 4_Compare

Figure 1 - Xenium prime (Ovarian cancer)
- 1_Cell_typing_Ovarian.ipynb
- 2_infercnv_Ovarian_moments_correct.ipynb
- 3_infercnv_ovarian_moments_ACROSS_winsize_n_read_depth.ipynb
- 4_spatial_stats_Ovarian_moments.ipynb

Figure 1 - Xenium prime (Lymph node)
- 1_Cell_typing_lymphnode.ipynb
- 2_infercnv_lymphnode_moments.ipynb
