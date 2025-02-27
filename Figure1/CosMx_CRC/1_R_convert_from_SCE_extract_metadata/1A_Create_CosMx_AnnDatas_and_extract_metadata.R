
# Load necessary libraries
library(dplyr)
library(HDF5Array)
library(zellkonverter)
library(stringr)




### Save the cluster annotations as metadata for each sample 

for (SID in c(
  '110',
  '120'
  '221', 
  '222', 
  '231', 
  '232', 
  '242'
  )){
  lv1 <- readRDS(paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/lv1-",SID,".rds"))
  cluster_metadata <- as.data.frame(lv1[["clust"]])
  write.csv(cluster_metadata, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",SID,"/cluster_metadata.csv"))
}


### Save the raw SCE as an AnnData object .h5ad file

raw_paths <- c(
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-110-20241210T191725Z-001/raw-110",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-120-20241210T191724Z-002/raw-120"
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-221-20241111T094332Z-001/raw-221",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-222-20241111T094333Z-001/raw-222",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-231-20241111T095538Z-001/raw-231",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-232-20241108T073407Z-001/raw-232",
  "~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-242-20241108T073409Z-001/raw-242"
)

for (path_obj in raw_paths){
  SID <- str_sub(path_obj, -3)
  sce_from_h5 <- HDF5Array::loadHDF5SummarizedExperiment(path_obj)
  h5ad_path = paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",SID,"/adata_",SID,"_raw.h5ad")
  writeH5AD(sce_from_h5, 
            file=h5ad_path)
}


### Save ROI healthy cell annotations

ROI <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/roi.rds")

#, '222', '231', '232', '242'
for (SID in c(
  '110',
  '120'
  '221', 
  '222', 
  '231', 
  '232', 
  '242'
  )){
  write.csv(ROI[[SID]],
            file = paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",SID,"/ROI_",SID,"_normal.csv"),
            row.names = TRUE
  )
}








### Extracting the ROI (REF, TVA, CRC) for each sample

ROI <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/cells_by_sid_roi.rds")

# Initialize an empty list to store data frames for each sample
sample_data_frames <- list()

# Loop through each sample in the RDS data
for (sample in names(ROI)) {
  sample_rois <- ROI[[sample]]
  
  # Initialize a data frame to hold cell IDs and their associated ROIs
  sample_df <- data.frame()
  
  # Loop through each ROI for the sample
  for (roi in names(sample_rois)) {
    # Extract the cell IDs for the current ROI
    cell_ids <- sample_rois[[roi]]
    
    # Modify the ROI name to keep only the part after the first underscore
    roi<- sub("^[^_]*_", "", roi)
    
    # Create a data frame with cell IDs as row names and the ROI as a column
    roi_df <- data.frame(ROI = rep(roi, length(cell_ids)), row.names = cell_ids)
    
    # Append to the sample's data frame
    sample_df <- rbind(sample_df, roi_df)
  }
  
  # Store the sample's data frame in the list
  sample_data_frames[[sample]] <- sample_df
  
  # Write the sample's data frame to a CSV file
  write.csv(sample_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",sample, "/",sample,"_ROI.csv"), row.names = TRUE)
}


## Covert the features to csv
features_all <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_all.rds")
features_all_df <- as.data.frame(features_all)
write.csv(features_all_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_all.csv"))


features_epi <- readRDS("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_epi.rds")
features_epi_df <- as.data.frame(features_epi)
write.csv(features_epi_df, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/features_epi.csv"))



