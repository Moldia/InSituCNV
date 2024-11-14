library(SingleCellExperiment)
library(zellkonverter)

######### Epithelial

# Sample snPATHO-seq B-1970164 epi
sce_snPATHOseq_B197 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/lab_lv2-B-1970164-01-04_2,epi.rds")

writeH5AD(
  sce_snPATHOseq_B197,               # The SCE object
  file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/adata_lab_lv2-B-1970164-01-04_2,epi.h5ad",
  X_name = NULL,                     # Use default primary matrix
  skip_assays = FALSE,               # Don't skip assays
  compression = "gzip",              # Use gzip compression (or choose "none" or "lzf")
  verbose = TRUE,                    # Print progress messages
  assays = TRUE,                     # Convert all assay data
  colData = TRUE,                    # Convert all colData (cell metadata)
  rowData = TRUE,                    # Convert all rowData (gene metadata)
  varm = TRUE,                       # Convert varm
  reducedDims = TRUE,                # Convert reduced dimensions (PCA, UMAP, etc.)
  metadata = TRUE,                   # Convert all metadata
  colPairs = TRUE,                   # Convert all colPairs
  rowPairs = TRUE                    # Convert all rowPairs
)

# Sample snPATHO-seq B-2080151 epi
snPATHOseq_B208 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/lab_lv2-B-2080151-01-18_2,epi.rds")

writeH5AD(
  snPATHOseq_B208,                   # The SCE object
  file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/adata_lab_lv2-B-2080151-01-18_2,epi.h5ad",
  X_name = NULL,                     # Use default primary matrix
  skip_assays = FALSE,               # Don't skip assays
  compression = "gzip",              # Use gzip compression (or choose "none" or "lzf")
  verbose = TRUE,                    # Print progress messages
  assays = TRUE,                     # Convert all assay data
  colData = TRUE,                    # Convert all colData (cell metadata)
  rowData = TRUE,                    # Convert all rowData (gene metadata)
  varm = TRUE,                       # Convert varm
  reducedDims = TRUE,                # Convert reduced dimensions (PCA, UMAP, etc.)
  metadata = TRUE,                   # Convert all metadata
  colPairs = TRUE,                   # Convert all colPairs
  rowPairs = TRUE                    # Convert all rowPairs
)

##### Immune

# Sample snPATHO-seq B-1970164 immune
sce_snPATHOseq_B197 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/lab_lv2-B-1970164-01-04_2,imm.rds")

writeH5AD(
  sce_snPATHOseq_B197,               # The SCE object
  file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/adata_lab_lv2-B-1970164-01-04_2,imm.h5ad",
  X_name = NULL,                     # Use default primary matrix
  skip_assays = FALSE,               # Don't skip assays
  compression = "gzip",              # Use gzip compression (or choose "none" or "lzf")
  verbose = TRUE,                    # Print progress messages
  assays = TRUE,                     # Convert all assay data
  colData = TRUE,                    # Convert all colData (cell metadata)
  rowData = TRUE,                    # Convert all rowData (gene metadata)
  varm = TRUE,                       # Convert varm
  reducedDims = TRUE,                # Convert reduced dimensions (PCA, UMAP, etc.)
  metadata = TRUE,                   # Convert all metadata
  colPairs = TRUE,                   # Convert all colPairs
  rowPairs = TRUE                    # Convert all rowPairs
)


# Sample snPATHO-seq B-2080151 immune
snPATHOseq_B208 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/lab_lv2-B-2080151-01-18_2,imm.rds")

writeH5AD(
  snPATHOseq_B208,                   # The SCE object
  file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/adata_lab_lv2-B-2080151-01-18_2,imm.h5ad",
  X_name = NULL,                     # Use default primary matrix
  skip_assays = FALSE,               # Don't skip assays
  compression = "gzip",              # Use gzip compression (or choose "none" or "lzf")
  verbose = TRUE,                    # Print progress messages
  assays = TRUE,                     # Convert all assay data
  colData = TRUE,                    # Convert all colData (cell metadata)
  rowData = TRUE,                    # Convert all rowData (gene metadata)
  varm = TRUE,                       # Convert varm
  reducedDims = TRUE,                # Convert reduced dimensions (PCA, UMAP, etc.)
  metadata = TRUE,                   # Convert all metadata
  colPairs = TRUE,                   # Convert all colPairs
  rowPairs = TRUE                    # Convert all rowPairs
)



##### Stromal

# Sample snPATHO-seq B-1970164 stromal
sce_snPATHOseq_B197 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/lab_lv2-B-1970164-01-04_2,str.rds")

writeH5AD(
  sce_snPATHOseq_B197,               # The SCE object
  file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/adata_lab_lv2-B-1970164-01-04_2,str.h5ad",
  X_name = NULL,                     # Use default primary matrix
  skip_assays = FALSE,               # Don't skip assays
  compression = "gzip",              # Use gzip compression (or choose "none" or "lzf")
  verbose = TRUE,                    # Print progress messages
  assays = TRUE,                     # Convert all assay data
  colData = TRUE,                    # Convert all colData (cell metadata)
  rowData = TRUE,                    # Convert all rowData (gene metadata)
  varm = TRUE,                       # Convert varm
  reducedDims = TRUE,                # Convert reduced dimensions (PCA, UMAP, etc.)
  metadata = TRUE,                   # Convert all metadata
  colPairs = TRUE,                   # Convert all colPairs
  rowPairs = TRUE                    # Convert all rowPairs
)


# Sample snPATHO-seq B-2080151 stromal
snPATHOseq_B197 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/lab_lv2-B-2080151-01-18_2,str.rds")

writeH5AD(
  snPATHOseq_B197,                   # The SCE object
  file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/snPATHO-seq/adata_lab_lv2-B-2080151-01-18_2,str.h5ad",
  X_name = NULL,                     # Use default primary matrix
  skip_assays = FALSE,               # Don't skip assays
  compression = "gzip",              # Use gzip compression (or choose "none" or "lzf")
  verbose = TRUE,                    # Print progress messages
  assays = TRUE,                     # Convert all assay data
  colData = TRUE,                    # Convert all colData (cell metadata)
  rowData = TRUE,                    # Convert all rowData (gene metadata)
  varm = TRUE,                       # Convert varm
  reducedDims = TRUE,                # Convert reduced dimensions (PCA, UMAP, etc.)
  metadata = TRUE,                   # Convert all metadata
  colPairs = TRUE,                   # Convert all colPairs
  rowPairs = TRUE                    # Convert all rowPairs
)

