library(dplyr)


# Sample WTx-CosMX B-1970164

WTcosMX_B197 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/WTx-CosMx/add-B-1970164-01-04_2.rds")

cell_metadata <- as.data.frame(colData(WTcosMX_B197))

selected_columns <- c("cell_id", "ist", "sub")  # Replace with your actual column names
selected_metadata <- cell_metadata %>% 
  select(all_of(selected_columns))

write.csv(selected_metadata, file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/WTx-CosMx/metadata_add-B-1970164-01-04_2.csv", row.names = TRUE)


raw_221  <- readRDS("/home/augusta/storage3/insituCNV/data/WTx-CosMx_TVA/round2/raw-221-20241111T094332Z-001/raw-221/se.rds")

# Sample WTx-CosMX B-2080151

WTcosMX_B208 <- readRDS("~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/WTx-CosMx/add-B-2080151-01-18_2.rds")

cell_metadata <- as.data.frame(colData(WTcosMX_B208))

selected_columns <- c("cell_id", "ist", "sub")  # Replace with your actual column names
selected_metadata <- cell_metadata %>% 
  select(all_of(selected_columns))

write.csv(selected_metadata, file = "~/storage3_mount/insituCNV/data/WTx-CosMx_TVA/WTx-CosMx/metadata_add-B-2080151-01-18_2.csv", row.names = TRUE)
