#################################################
# RNA-seq Data Preprocessing Script
# Project: E. coli RmpA Regulation Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)

# Create output directories if they don't exist
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/qc", recursive = TRUE, showWarnings = FALSE)

###################
# Data Import 
###################
# Import count data
counts_data <- read.delim("data/raw/GSE286114_GPL21433_series_matrix.txt.gz", 
                         skip = 0)

# Create sample metadata
sample_metadata <- data.frame(
    sample_id = c("GSM8717926", "GSM8717927", "GSM8717928",
                 "GSM8717929", "GSM8717930", "GSM8717931"),
    condition = factor(c(rep("control", 3), rep("rmpA_overexpression", 3))),
    replicate = rep(1:3, 2)
)

# Save raw metadata
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

###################
# Create DESeq2 Object
###################
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = sample_metadata,
    design = ~ condition
)

# Save raw counts
saveRDS(counts_data, "data/processed/raw_counts.rds")

###################
# QC Plots
###################
pdf("figures/qc/raw_data_qc.pdf")

# Library size plot
barplot(colSums(counts(dds)), 
        main="Library Sizes",
        names=colnames(counts(dds)),
        las=2)

# Genes detected per sample
expressed_genes <- colSums(counts(dds) > 0)
barplot(expressed_genes,
        main="Number of Expressed Genes",
        names=colnames(counts(dds)),
        las=2)

dev.off()

###################
# Normalization
###################
# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Save normalized data
saveRDS(dds, "data/processed/dds_object.rds")
saveRDS(normalized_counts, "data/processed/normalized_counts.rds")

# Save processing info
sessionInfo()
