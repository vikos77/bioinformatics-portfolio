# Script to download and prepare E. coli RNA-seq data
library(GEOquery)
library(tidyverse)

# Download expression data and metadata from GEO
gse <- getGEO("GSE115313", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

# Save metadata
write_csv(metadata, "data/raw/metadata.csv")

# Download supplementary files
getGEOSuppFiles("GSE115313", baseDir = "data/raw")
