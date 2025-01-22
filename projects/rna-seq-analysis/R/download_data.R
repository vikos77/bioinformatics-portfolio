# Script to download and prepare bacterial RNA-seq data for RmpA regulation study
library(GEOquery)
library(tidyverse)

# Download expression data and metadata from GEO
gse <- getGEO("GSE286114", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

# Save metadata
write_csv(metadata, "data/raw/metadata.csv")

# Download supplementary files
getGEOSuppFiles("GSE286114", baseDir = "data/raw")
