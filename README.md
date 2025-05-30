# InSituCNV
This is a public repository to reproduce the analysis in the manusctipt ***In Situ* Inference of Copy Number Variations in Image-Based Spatial Transcriptomics** by Jensen et al., *unpublished*.

link: 

***
## Data avaliability
Raw data of the CosMX WT experiments in CRC samples included in this work can be downloaded from **ADD LINK**. The 10X Xenium prime datasets used in this study are publically avaliable at the following link: **ADD LINK**. The 10X scRNA-seq dataset used for CNV simulation is avaliable at **ADD LINK**.

From Cowell, L.G., Pascual Reguant, A., Heyn, H., *et al.* (2025). *[Title of the preprint]*. bioRxiv.:
- CosMx WTx: **ADD LINK**
- snPATHO-seq: **ADD LINK**

From 10X publically avaliable datasets:
- Xenium prime (Ovarian): **ADD LINK**
- Xenium prime (lymph node): **ADD LINK**
- scRNAseq (lung organoid): **ADD LINK**

The data is avalible 

*** 
## Folder structure
These folders contains all R script and notebooks necessary to reproduce the analysis in this study. For simplicity, they are organized in the same order as they are presented in the manuscipt.

### Figure 1
- Figure 1 - CosMx and snPATHO-seq (Colorectal cancer)
  - 1_R_convert_from_SCE_extract_metadata
  - 2_Prepare_data
  - 3_inferCNV
  - 4_Compare
- Figure 1 - Xenium prime (Ovarian cancer)
  - 01_
  - 02_
  - 03_
  - 04_
- Figure 1 - Xenium prime (Lymph node)
  - 1_Cell_typing_lymphnode.ipynb
  - 2_infercnv_lymphnode_moments.ipynb

### Figure 2
- 00_Load_simulation_dataset
- 01_Simulate_CNVs
- 02_Apply_technical_variations
- 03_Run_inferCNV
- 04_Compare_results_using_metrics
- 05_Visualize_results

### insitucnv
Package containing the functions used in this study. 


### data (?)
Location of the input data files necessary to run the analysis. Intermediate data files have to be generated from the input files. 



***
## How to use:
Clone the repository:

```git clone https://github.com/Moldia/InSituCNV.git```

Navigate to the directory:

```cd InSituCNV```

Create and activate a conda environment including all necessary dependencies:

```conda env create --name insituCNV_env --file=insitucnv.yml```

```conda activate insituCNV_env```



