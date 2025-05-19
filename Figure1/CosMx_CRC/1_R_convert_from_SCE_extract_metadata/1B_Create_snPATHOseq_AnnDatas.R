library(SingleCellExperiment)
library(zellkonverter)
library(DropletUtils)   # For reading 10x Genomics data
library(stringr)


snPATHOseq <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/snPATHO-seq.rds")

duplicated_names <- colnames(colData(snPATHOseq))
duplicated_names[duplicated(duplicated_names)]

colnames(colData(snPATHOseq)) <- make.unique(colnames(colData(snPATHOseq)))



writeH5AD(
  snPATHOseq,               # The SCE object
  file = "/home/augusta/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/snPATHOseq.h5ad",
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




#For the raw snPATHO-seq data

raw_paths <- c(
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/110-20250209T170643Z-001/110",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/120-20250209T170642Z-001/120",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/210-20250209T170645Z-001/210",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/221-20250209T170647Z-001/221"
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/222-20250209T170637Z-001/222",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/230-20250209T170640Z-001/230",
  #"~/storage3//insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/242-20250209T170638Z-001/242"
)

for (path_obj in raw_paths){
  SID <- str_sub(path_obj, -3)
  snPATHOseq <- read10xCounts(path_obj)
  snPATHOseq$Sample <- SID
  
  writeH5AD(
    snPATHOseq,               # The SCE object
    file = paste0("/home/augusta/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/snPATHOseq_raw_",SID,".h5ad"),
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
}
