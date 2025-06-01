# Load necessary libraries
library(dplyr)           
library(HDF5Array)       
library(zellkonverter)  
library(stringr)       


###------------------------------------------------------
### Save the raw SCE as an AnnData object (.h5ad file)

# Define paths to raw HDF5-based SummarizedExperiment objects
raw_paths <- c(
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-110-20241210T191725Z-001/raw-110",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-120-20241210T191724Z-002/raw-120"
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-221-20241111T094332Z-001/raw-221",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-222-20241111T094333Z-001/raw-222",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-231-20241111T095538Z-001/raw-231",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-232-20241108T073407Z-001/raw-232",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-242-20241108T073409Z-001/raw-242",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-210-20250331T141016Z-001/raw-210",
  #"~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-241-20250408T192509Z-001/raw-241"
)

# Iterate through each path object
for (path_obj in raw_paths) {
  # Extract sample ID (last 3 characters of the path)
  SID <- str_sub(path_obj, -3)
  
  # Define output path for the .h5ad file
  h5ad_path = paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/", SID, "/adata_", SID, "_raw.h5ad")
  
  # Load the HDF5-based SummarizedExperiment and save it as a .h5ad file
  writeH5AD(HDF5Array::loadHDF5SummarizedExperiment(path_obj), file = h5ad_path)
}


###------------------------------------------------------
### Save the filtered cells' cluster annotations as metadata for each sample

# Loop through specified sample ID
for (SID in c('221')) {
  # Load level 1 clustering RDS file for the sample
  lv1 <- readRDS(paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/lv1-", SID, ".rds"))
  
  # Extract clustering information and convert to data frame
  cluster_metadata <- as.data.frame(lv1[["clust"]])
  
  # Save the cluster metadata as a CSV file within the sample directory
  write.csv(cluster_metadata, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/", SID, "/cluster_metadata.csv"))
}

###------------------------------------------------------
### Extract and save per-sample healthy ROI (REF) annotations

ROI <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/roi.rds")

#, '222', '231', '232', '242'
for (SID in c(
  '210'
  #'110',
  #'120'
  #'221', 
  #'222', 
  #'231', 
  #'232', 
  #'242'
  )){
  write.csv(ROI[[SID]],
            file = paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",SID,"/ROI_",SID,"_normal.csv"),
            row.names = TRUE
  )
}


###------------------------------------------------------
### Extract and save per-sample ROI (REF, TVA1-3, CRC) annotations

# Load RDS file mapping samples to ROI-associated cell IDs
ROI <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/cells_by_sid_roi.rds")

# Initialize a list to store output data frames for each sample
sample_data_frames <- list()

# Loop over all samples in the ROI object
for (sample in names(ROI)) { 
  sample_rois <- ROI[[sample]]
  
  # Prepare an empty data frame for this sample
  sample_df <- data.frame()
  
  # Iterate over each ROI within the sample
  for (roi in names(sample_rois)) {
    # Get cell IDs associated with the current ROI
    cell_ids <- sample_rois[[roi]]
    
    # Simplify ROI label by removing text before the first underscore
    roi <- sub("^[^_]*_", "", roi)
    
    # Create a data frame with cell IDs and corresponding ROI label
    roi_df <- data.frame(ROI = rep(roi, length(cell_ids)), row.names = cell_ids)
    
    # Append to the cumulative sample data frame
    sample_df <- rbind(sample_df, roi_df)
  }
  
  # Save the sample's ROI annotation to list and CSV
  sample_data_frames[[sample]] <- sample_df
  write.csv(sample_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/", sample, "/", sample, "_ROI.csv"), row.names = TRUE)
}


###------------------------------------------------------
### Convert the features list to a CSV files

features_all <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_all.rds")
features_all_df <- as.data.frame(features_all)
write.csv(features_all_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_all.csv"))

features_epi <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_epi.rds")
features_epi_df <- as.data.frame(features_epi)
write.csv(features_epi_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_epi.csv"))

epi_gs <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/epi_gs.rds")
epi_gs_df <- as.data.frame(epi_gs)
write.csv(epi_gs_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/epi_gs.csv"))

