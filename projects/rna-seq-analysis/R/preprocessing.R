#################################################
# RNA-seq Data Preprocessing Script
# Project: E. coli RmpA Regulation Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(GEOquery)
library(readxl)
library(R.utils)

# Set relative paths for output directories (relative to script location)
output_data_dir <- "../data/processed"
output_qc_results_dir <- "../results/qc"
output_qc_figures_dir <- "../figures/qc"
output_raw_data_dir <- "../data/raw"

# Create output directories if they don't exist
dir.create(output_raw_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_qc_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_qc_figures_dir, recursive = TRUE, showWarnings = FALSE)


###################
# Data Download
###################

# Download the dataset form GEO
gse <- getGEO("GSE286114", GSEMatrix = TRUE, destdir = output_raw_data_dir) # Save raw data files to data/raw folder
metadata <- pData(phenoData(gse[[1]]))

# Get supplementary files
getGEOSuppFiles("GSE286114", baseDir = output_raw_data_dir) # Save supplementary files to data/raw folder

###################
# Data Import
###################

# Create temporary file and decompress
temp_file <- tempfile(fileext = ".xls")

raw_fpkm_file_gz <- file.path(output_raw_data_dir, "GSE286114", "GSE286114_gene_fpkm_RmpA_Escherichia.xls.gz")

if (file.exists(raw_fpkm_file_gz)) {
  R.utils::gunzip(raw_fpkm_file_gz,
                  destname = temp_file,
                  remove = FALSE)
} else {
  stop("File not found: Ensure GEO download was successful. Expected file: ", raw_fpkm_file_gz)
}


# Read the Excel file
counts_data <- read_excel(temp_file)


# Examine the data structure
dim(counts_data)  # Check dimensions
head(counts_data) # Look at first few rows
colnames(counts_data) # Check column names

# Clean up column names and extract expression data
expression_data <- counts_data %>%
  select(gene_id, PLB1K1, PLB1K2, PLB1K3, `PLB1K-rmpA1`, `PLB1K-rmpA2`, `PLB1K-rmpA3`) %>%
  mutate(across(-gene_id, as.numeric)) %>%  # Convert all columns except gene_id to numeric
  column_to_rownames("gene_id")

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = c("PLB1K1", "PLB1K2", "PLB1K3", "PLB1K-rmpA1", "PLB1K-rmpA2", "PLB1K-rmpA3"),
  condition = factor(c(rep("control", 3), rep("rmpA_overexpression", 3))),
  replicate = rep(1:3, 2),
  row.names = c("PLB1K1", "PLB1K2", "PLB1K3", "PLB1K-rmpA1", "PLB1K-rmpA2", "PLB1K-rmpA3")
)

# Save sample metadata to processed data directory
saveRDS(sample_metadata, file.path(output_data_dir, "sample_metadata.rds"))

#Save gene information separately
gene_info <- counts_data %>%
  select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand, gene_length, gene_biotype, gene_description)

write_csv(gene_info, file.path(output_data_dir, "gene_info.csv"))

###################
# Create DESeq2 Object
###################

#Convert FPKM to counts (approximately)
# Round to nearest integer as DESeq2 requires count data
counts_matrix <- round(expression_data*1000)

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_metadata,
  design = ~ condition
)

#Look at the distribution of counts
total_counts_per_gene <- rowSums(counts(dds))
summary(total_counts_per_gene)

# Create a histogram of the count distribution and save to figures/qc
png(file.path(output_qc_figures_dir, "gene_count_distribution_before_filtering.png"))
hist(log10(total_counts_per_gene + 1),
     breaks = 50,
     main = "Distribution of Total Counts per Gene (Before Filtering)",
     xlab = "log10(Total Counts + 1)")
dev.off()


#Apply filtering
#Filtering is applied to remove noise from our data which might affect the statistical analysis.

# Data-driven filtering
min_count <- 10  # minimum counts per sample
min_samples <- 3 # minimum number of samples
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds_filtered <- dds[keep,]

# Create a histogram of the count distribution and save to figures/qc
png(file.path(output_qc_figures_dir, "gene_count_distribution_after_filtering.png"))
hist(log10(total_counts_per_gene + 1),
     breaks = 50,
     main = "Distribution of Total Counts per Gene (After Filtering)",
     xlab = "log10(Total Counts + 1)")

# Save filtering statistics to results/qc
filtering_stats <- data.frame(
  total_genes = nrow(dds),
  genes_kept = sum(keep),
  genes_filtered = nrow(dds) - sum(keep)
)
print(filtering_stats)
write.csv(filtering_stats, file.path(output_qc_results_dir, "filtering_statistics.csv"))

# Visualize counts per sample before and after filtering - save to figures/qc
png(file.path(output_qc_figures_dir, "count_boxplots_filtering_effect.png"))
par(mfrow=c(1,2))

# Boxplot of raw counts (before filtering)
boxplot(log2(counts(dds) + 1),
        main="Raw counts",
        ylab="log2(counts)",
        las=2,
        )

# Boxplot of counts after filtering
boxplot(log2(counts(dds_filtered) + 1),
        main="After Filtering",
        ylab="log2(counts + 1)",
        las=2,
        )
dev.off()

# Save filtered DESeq object to processed data directory
saveRDS(dds_filtered, file.path(output_data_dir, "dds_filtered_object.rds"))

###################
# QC Plots
###################

# Create QC plots for filtered data and save to figures/qc
pdf(file.path(output_qc_figures_dir, "filtered_data_qc.pdf"))

# Sample correlation heatmap
vsd <- vst(dds_filtered, blind=TRUE)
sampleDists <- dist(t(assay(vsd)))

#Heatmap showing sample-to-sample distances based on gene expression. Used to check for overall sample relationships and potential outliers.
pheatmap(as.matrix(sampleDists),
         main="Sample Distance Matrix (Filtered Data)",
         )

# Enhanced PCA plot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Principal Component Analysis (PCA) plot of filtered and variance-stabilized data. Used to assess sample clustering by experimental condition."
print(
  ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot of Filtered Data") +
    theme_minimal
)

dev.off()


###################
# Normalization
###################
# Now proceed with normalization
dds_filtered <- estimateSizeFactors(dds_filtered)
normalized_counts <- counts(dds_filtered, normalized=TRUE)

# Save normalized data to processed data directory
saveRDS(normalized_counts, file.path(output_data_dir, "normalized_counts.rds"))

# Create boxplots to visualize normalization effect
png(file.path(output_qc_figures_dir, "count_boxplots_normalization_effect.png"))
par(mfrow=c(1,2))


#Boxplot of log2-transformed raw counts (filtered data) before normalization. 
#Shows potential differences in count distributions across samples.
boxplot(log2(counts(dds_filtered)),
        main="Raw counts",
        ylab="log2(counts)",
        las=2,
        )

# Boxplot of log2-transformed normalized counts.
#Demonstrates the effect of normalization in making count distributions more comparable across samples.
boxplot(log2(normalized_counts),
        main="Normalized counts",
        ylab="log2(counts)",
        las=2,
        )

dev.off()

# Save normalized counts to processed data directory (redundant, already saved above, but kept for clarity)
saveRDS(normalized_counts, file.path(output_data_dir, "normalized_counts.rds"))

# Save processing info (session info) to results/qc
capture.output(sessionInfo(), file = file.path(output_qc_results_dir, "session_info.txt"))