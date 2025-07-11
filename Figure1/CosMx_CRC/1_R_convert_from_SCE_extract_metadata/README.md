## R Dependencies

This project uses R version **4.5.1** with the following core packages:

- `dplyr`
- `stringr`
- `HDF5Array` 
- `zellkonverter` 

To install the required R packages, run the following in an R session:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("HDF5Array", "zellkonverter", "dplyr", "stringr"))
