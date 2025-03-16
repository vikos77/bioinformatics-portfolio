#################################################
# RNA-seq Exploratory Analysis Script
# Project: Fibrosis Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
})

# Logging function
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting exploratory analysis")

# Load the DESeq2 object from preprocessing
dds <- readRDS("results/dds_object.rds")

###################
# Data Normalization
###################

log_message("Performing normalization")

# Calculate size factors to account for library size differences
dds <- estimateSizeFactors(dds)

# Display the size factors for each sample
sizeFactors(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

# Compare raw vs normalized counts for an example gene
gene_id <- "ENSMUSG00000051951"  
raw_counts <- counts(dds)[gene_id,]
norm_counts <- normalized_counts[gene_id,]
comparison <- data.frame(
  Sample = names(raw_counts),
  Raw = raw_counts,
  Normalized = round(norm_counts, 2)
)
print(comparison)

# Save normalized counts for later use
saveRDS(normalized_counts, "results/normalized_counts.rds")

###################
# Variance Stabilizing Transformation
###################

log_message("Performing variance stabilizing transformation")

# Apply variance stabilizing transformation for visualization
# blind=TRUE prevents the design formula from influencing the transformation
vsd <- vst(dds, blind = TRUE)

###################
# Sample Correlation Analysis
###################

log_message("Creating sample correlation heatmap")

# Extract the transformed count matrix
vsd_mat <- assay(vsd) 

# Calculate correlation matrix between samples
vsd_cor <- cor(vsd_mat) 

# Create directory for figures if it doesn't exist
dir.create("results/figures/exploratory", recursive = TRUE, showWarnings = FALSE)

# Save correlation matrix as CSV
write.csv(vsd_cor, "results/sample_correlations.csv")

# Basic correlation heatmap
png("results/figures/exploratory/correlation_heatmap_basic.png", width=480, height=480)
pheatmap(vsd_cor)
dev.off()

# Enhanced correlation heatmap with condition annotation
png("results/figures/exploratory/correlation_heatmap_annotated.png", width=480, height=480)
pheatmap(vsd_cor, 
         annotation = data.frame(Condition = dds$condition, 
                                 row.names = colnames(dds)),
         annotation_colors = list(Condition = c(fibrosis = "red", normal = "blue")),
         main = "Sample Correlation by Condition")
dev.off()

###################
# Principal Component Analysis
###################

log_message("Performing PCA analysis")

# Basic PCA plot
png("results/figures/exploratory/pca_basic.png", width=480, height=480)
plotPCA(vsd, intgroup="condition")
dev.off()

# Extract PCA data for customization
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Save PCA data for reference
saveRDS(pca_data, "results/pca_data.rds")

# Enhanced PCA plot with better styling
png("results/figures/exploratory/pca_enhanced.png", width=480, height=480)
ggplot(pca_data, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principal Component Analysis") +
  theme_bw() +
  scale_color_manual(values = c(fibrosis = "red", normal = "blue")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
dev.off()

log_message("Exploratory analysis complete")
