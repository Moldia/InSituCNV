
# Load necessary libraries
library(dplyr)
library(HDF5Array)
library(zellkonverter)
library(stringr)




### Save the cluster annotations as metadata for each sample 

for (SID in c('221', '222', '231', '232', '242')){
  lv1 <- readRDS(paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/lv1-",SID,".rds"))
  cluster_metadata <- as.data.frame(lv1[["clust"]])
  write.csv(cluster_metadata, paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",SID,"/cluster_metadata.csv"))
}




### Save the raw SCE as an AnnData object .h5ad file

raw_paths <- c(
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
for (SID in c('221', '222', '231', '232', '242')){
  write.csv(ROI[[SID]],
            file = paste0("~/storage3/insituCNV/data/WTx-CosMx_TVA/round2/",SID,"/ROI_",SID,"_normal.csv"),
            row.names = TRUE
  )
}
