# InSituCNV
This is a public repository to reproduce the analysis in the manusctipt ***In Situ* Inference of Copy Number Variations in Image-Based Spatial Transcriptomics** by Jensen et al., https://doi.org/10.1101/2025.07.02.662761.

*** 
## Code and Data Avaliability
These folders contains all R script and notebooks necessary to reproduce the analysis in this study. For simplicity, they are organized in the same order as they are presented in the manuscipt. 

### Figure 1
#### CosMx_CRC
Includes the analysis performed on the CosMx and snPATHO-seq datasets from colorectal tumor sample 221, avaliable in Crowell et al., "*Tracing colorectal malignancy transformation from cell to tissue scale*", bioRxiv (2025) https://doi.org/10.1101/2025.06.23.660674.
#### Xenium_prime_ovarian
Contains analysis of the publicly available high-grade ovarian cancer Xenium Prime dataset from 10x Genomics: https://www.10xgenomics.com/datasets/xenium-prime-fresh-frozen-human-ovary.
#### Xenium_prime_lymph_node
Contains analysis of the publicly available lymph node Xenium Prime dataset from 10x Genomics: https://www.10xgenomics.com/datasets/preview-data-xenium-prime-gene-expression.

### Figure 2
This folder contains the code used to assess the technical limitations of CNV inference, as presented in Figure 2 of the manuscript "In Situ Inference of Copy Number Variations in Image-Based Spatial Transcriptomics". The 10x scRNA-seq dataset used for CNV simulation is available at: https://cellxgene.cziscience.com/collections/e9cf4e8d-05ed-4d95-b550-af35ca219390.



***
## How to use:
Clone the repository:

```git clone https://github.com/Moldia/InSituCNV.git```

Navigate to the directory:

```cd InSituCNV```

Create and activate a conda environment including all necessary dependencies:

```conda env create --name insituCNV_env --file=insitucnv.yml```

```conda activate insituCNV_env```



