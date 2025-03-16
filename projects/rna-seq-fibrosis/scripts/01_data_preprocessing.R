#################################################
# RNA-seq Data Preprocessing Script
# Project: Fibrosis Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

# Create a function to log the progress
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting fibrosis RNA-seq data preprocessing")

# Read the count data file
log_message("Reading raw count data")
wt_rawcounts <- read.csv("data/fibrosis_smoc2_rawcounts.csv")

# Check the structure of the data
str(wt_rawcounts)

# Convert the gene name column to row names 
wt_rawcounts <- wt_rawcounts %>% 
  column_to_rownames("X")

log_message("Preparing experimental metadata")

# Create metadata for the DESeq2 analysis
# Define sample genotypes
genotype <- c("wt", "wt", "wt", "wt", "wt", "wt", "wt")

# Define experimental conditions for each sample
condition <- c("fibrosis", "fibrosis", "normal", "normal", "fibrosis", "normal", "fibrosis")

# Combine into a dataframe for the analysis
wt_dataframe <- data.frame(genotype, condition)

# Add rownames to the metadata, ensuring they match column names in the raw count data
rownames(wt_dataframe) <- colnames(wt_rawcounts)

log_message("Creating DESeq2 object")

# Create the DESeq2 object by providing:
# 1. countData: Raw counts matrix with genes as rows and samples as columns
# 2. colData: Sample metadata with experimental variables
# 3. design: Formula specifying the model for testing (using condition as the factor)
dds <- DESeqDataSetFromMatrix(
  countData = wt_rawcounts,
  colData = wt_dataframe,
  design = ~ condition
)

# Verify conditions are correctly specified
levels(dds$condition)

# Summarize the DESeq2 object
summary(dds)

# Save the DESeq2 object for subsequent analysis
saveRDS(dds, "results/dds_object.rds")

log_message("Preprocessing complete. DESeq2 object saved.")
