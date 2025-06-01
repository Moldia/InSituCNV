# Load required libraries
library(SingleCellExperiment)  
library(zellkonverter)        
library(DropletUtils)       
library(stringr)               

###------------------------------------------------------
### Load and process the annotated snPATHO-seq SingleCellExperiment object

# Read the processed SCE object containing snPATHO-seq data
snPATHOseq <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/snPATHO-seq.rds")

# Check for duplicated column names in the colData (cell metadata)
duplicated_names <- colnames(colData(snPATHOseq))
duplicated_names[duplicated(duplicated_names)]  # Print duplicated column names

# Make all colData column names unique to avoid issues with export
colnames(colData(snPATHOseq)) <- make.unique(colnames(colData(snPATHOseq)))

# Export the SCE object to AnnData format (.h5ad) for use in Python/Scanpy
writeH5AD(
  snPATHOseq,               # The SingleCellExperiment object
  file = "/home/augusta/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/snPATHOseq.h5ad",
  X_name = NULL,            # Use the default assay as the main matrix
  skip_assays = FALSE,      # Include all assay data in export
  compression = "gzip",     # Apply gzip compression to reduce file size
  verbose = TRUE,           # Print progress to console
  assays = TRUE,            # Export all assay matrices
  colData = TRUE,           # Export all column metadata (cell-level)
  rowData = TRUE,           # Export all row metadata (gene-level)
  varm = TRUE,              # Export variable-level metadata if present
  reducedDims = TRUE,       # Export reduced dimensional embeddings (e.g., PCA, UMAP)
  metadata = TRUE,          # Export additional metadata stored in the SCE object
  colPairs = TRUE,          # Export column-based pairwise relationships if present
  rowPairs = TRUE           # Export row-based pairwise relationships if present
)

###------------------------------------------------------
### Convert raw snPATHO-seq data from 10x Genomics format to AnnData (.h5ad)

# Define paths to raw 10x Genomics output folders for each sample
raw_paths <- c(
  # Paths for other samples are commented out for now
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/110-20250209T170643Z-001/110",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/120-20250209T170642Z-001/120",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/210-20250209T170645Z-001/210",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/221-20250209T170647Z-001/221"
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/222-20250209T170637Z-001/222",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/230-20250209T170640Z-001/230",
  #"~/storage3//insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/242-20250209T170638Z-001/242"
)

# Loop through each raw data path and convert to .h5ad
for (path_obj in raw_paths){
  # Extract the 3-digit sample ID from the path
  SID <- str_sub(path_obj, -3)
  
  # Read the raw 10x Genomics count matrix as a SingleCellExperiment object
  snPATHOseq <- read10xCounts(path_obj)
  
  # Add sample ID to the cell metadata
  snPATHOseq$Sample <- SID
  
  # Export the raw snPATHO-seq object as a .h5ad file
  writeH5AD(
    snPATHOseq,
    file = paste0("/home/augusta/storage3/insituCNV/data/WTx-CosMx_TVA/round2/snPATHOseq/snPATHOseq_raw_", SID, ".h5ad"),
    X_name = NULL,
    skip_assays = FALSE,
    compression = "gzip",
    verbose = TRUE,
    assays = TRUE,
    colData = TRUE,
    rowData = TRUE,
    varm = TRUE,
    reducedDims = TRUE,
    metadata = TRUE,
    colPairs = TRUE,
    rowPairs = TRUE
  )
}
