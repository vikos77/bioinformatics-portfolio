################################################
# RNA-seq Data Preprocessing Script
# Project: E. coli RmpA Regulation Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(biomaRt)
library(ggplot2)

# Set working directory
setwd("/media/vicky/OS/Users/kakar/bioinformatics-portfolio/projects/rna-seq-analysis")

###################
# Data Import 
###################

# Read count data
counts_data <- read.delim("data/raw/GSE286114_GPL21433_series_matrix.txt.gz", 
                         skip = 0)  # Adjust skip based on file format

# Create sample metadata
sample_metadata <- data.frame(
    sample_id = colnames(counts_data),
    condition = factor(c(rep("control", 3), rep("rmpA_overexpression", 3))),
    replicate = rep(1:3, 2)
)

###################
# Quality Control
###################

# 1. Basic count statistics
count_stats <- data.frame(
    total_counts = colSums(counts_data),
    detected_genes = colSums(counts_data > 0)
)

# Save QC metrics
write_csv(count_stats, "results/qc_metrics.csv")

# 2. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = sample_metadata,
    design = ~ condition
)

# 3. Pre-filtering low count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 4. Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Save normalized counts
write.csv(normalized_counts, "results/normalized_counts.csv")

###################
# QC Visualizations
###################

# 1. Count distribution plots
pdf("figures/count_distribution.pdf")
par(mfrow=c(2,1))
boxplot(log2(counts_data + 1), main="Raw Counts Distribution")
boxplot(log2(normalized_counts + 1), main="Normalized Counts Distribution")
dev.off()

# 2. Sample correlation heatmap
vsd <- vst(dds, blind=TRUE)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf("figures/sample_correlation_heatmap.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main="Sample Distance Matrix")
dev.off()

# 3. PCA plot
pdf("figures/pca_plot.pdf")
plotPCA(vsd, intgroup="condition")
dev.off()

# Save preprocessed data for downstream analysis
saveRDS(dds, "data/processed/dds_object.rds")
