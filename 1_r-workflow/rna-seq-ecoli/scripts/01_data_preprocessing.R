#################################################
# RNA-seq Data Preprocessing Script
# Project: RmpA Regulation Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(GEOquery)
library(readxl)
library(R.utils)
library(here)

# Create a function to log the progress
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting RmpA RNA-seq data preprocessing")

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/qc", recursive = TRUE, showWarnings = FALSE)

###################
# Data Download and Import
###################

log_message("Downloading data from GEO (GSE286114)")

# Check if data already exists to avoid redownloading
if (!file.exists("data/raw/GSE286114_family.soft.gz")) {
  tryCatch({
    gse <- getGEO("GSE286114", GSEMatrix = TRUE, destdir = "data/raw")
    metadata <- pData(phenoData(gse[[1]]))
    
    # Get supplementary files
    getGEOSuppFiles("GSE286114", baseDir = "data/raw")
    log_message("Data download successful")
  }, error = function(e) {
    log_message(paste("ERROR: Data download failed:", e$message))
    stop("Data download failed. Please check your internet connection or try manually downloading from GEO.")
  })
} else {
  log_message("GEO data already exists, using cached version")
  gse <- getGEO("GSE286114", GSEMatrix = TRUE, destdir = "data/raw")
  metadata <- pData(phenoData(gse[[1]]))
}

# Import expression data
log_message("Importing expression data")

temp_file <- tempfile(fileext = ".xls")
raw_fpkm_file_gz <- "data/raw/GSE286114/GSE286114_gene_fpkm_RmpA_Escherichia.xls.gz"

if (file.exists(raw_fpkm_file_gz)) {
  tryCatch({
    R.utils::gunzip(raw_fpkm_file_gz, destname = temp_file, remove = FALSE)
    counts_data <- read_excel(temp_file)
    log_message(paste("Successfully imported data with", nrow(counts_data), "genes"))
  }, error = function(e) {
    log_message(paste("ERROR: File import failed:", e$message))
    stop("File import failed. Please check file integrity.")
  })
} else {
  stop("File not found: ", raw_fpkm_file_gz)
}

# Clean up and prepare expression data
expression_data <- counts_data %>%
  select(gene_id, PLB1K1, PLB1K2, PLB1K3, `PLB1K-rmpA1`, `PLB1K-rmpA2`, `PLB1K-rmpA3`) %>%
  mutate(across(-gene_id, as.numeric)) %>%
  column_to_rownames("gene_id")

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = c("PLB1K1", "PLB1K2", "PLB1K3", "PLB1K-rmpA1", "PLB1K-rmpA2", "PLB1K-rmpA3"),
  condition = factor(c(rep("control", 3), rep("rmpA_overexpression", 3)), 
                     levels = c("control", "rmpA_overexpression")),
  replicate = rep(1:3, 2),
  batch = "batch1",
  row.names = c("PLB1K1", "PLB1K2", "PLB1K3", "PLB1K-rmpA1", "PLB1K-rmpA2", "PLB1K-rmpA3")
)

# Save sample metadata
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

# Save gene information
gene_info <- counts_data %>%
  select(gene_id, gene_name, gene_chr, gene_start, gene_end, 
         gene_strand, gene_length, gene_biotype, gene_description)

write_csv(gene_info, "data/processed/gene_info.csv")
log_message("Saved gene annotation information")

###################
# Create DESeq2 Object and Filtering
###################

log_message("Creating DESeq2 object")

# Convert FPKM to counts (approximately) and round to nearest integer
counts_matrix <- round(expression_data * 1000)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_metadata,
  design = ~ condition
)

# Look at count distribution before filtering
total_counts_per_gene <- rowSums(counts(dds))
log_message(paste("Count summary before filtering:",
                  "Min:", min(total_counts_per_gene),
                  "Median:", median(total_counts_per_gene),
                  "Max:", max(total_counts_per_gene)))

# Create count distribution plots
pdf("results/figures/qc/count_distributions.pdf", width=10, height=8)
par(mfrow=c(2,2))

# Plot raw count distribution
hist(log10(total_counts_per_gene + 1),
     breaks = 50,
     main = "Distribution of Total Counts per Gene\n(Before Filtering)",
     xlab = "log10(Total Counts + 1)",
     col = "steelblue",
     border = "white")

# Plot sample counts
boxplot(log10(counts(dds) + 1),
        main = "Log10 Counts per Sample\n(Before Filtering)",
        ylab = "log10(counts + 1)",
        las = 2,
        col = "lightblue",
        cex.axis = 0.8,
        names = gsub("PLB1K-", "R", colnames(dds)))

# Apply filtering
log_message("Applying count filtering")
min_count <- 10
min_samples <- 3
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds_filtered <- dds[keep,]

# Log filtering results
log_message(paste("Filtering results:",
                  "Total genes:", nrow(dds),
                  "Genes kept:", sum(keep),
                  "Genes filtered out:", nrow(dds) - sum(keep),
                  paste0("(", round((nrow(dds) - sum(keep))/nrow(dds) * 100, 2), "%)")))

# Plot filtered count distribution
total_counts_filtered <- rowSums(counts(dds_filtered))
hist(log10(total_counts_filtered + 1),
     breaks = 50,
     main = "Distribution of Total Counts per Gene\n(After Filtering)",
     xlab = "log10(Total Counts + 1)",
     col = "darkgreen",
     border = "white")

# Plot filtered sample counts
boxplot(log10(counts(dds_filtered) + 1),
        main = "Log10 Counts per Sample\n(After Filtering)",
        ylab = "log10(counts + 1)",
        las = 2,
        col = "lightgreen",
        cex.axis = 0.8,
        names = gsub("PLB1K-", "R", colnames(dds_filtered)))

dev.off()
log_message("Generated count distribution plots")

###################
# Normalization
###################

log_message("Normalizing count data")

# Perform normalization
dds_filtered <- estimateSizeFactors(dds_filtered)
normalized_counts <- counts(dds_filtered, normalized = TRUE)

# Create normalization effect plots
pdf("results/figures/qc/normalization_effect.pdf", width=10, height=6)
par(mfrow = c(1, 2))

# Boxplot of raw counts
boxplot(
  log2(counts(dds_filtered)),
  main = "Raw counts (filtered)",
  ylab = "log2(counts)",
  las = 2,
  col = "lightblue",
  cex.axis = 0.7,
  names = gsub("PLB1K-", "R", colnames(dds_filtered))
)

# Boxplot of normalized counts
boxplot(
  log2(normalized_counts),
  main = "Normalized counts",
  ylab = "log2(counts)",
  las = 2,
  col = "lightgreen",
  cex.axis = 0.7,
  names = gsub("PLB1K-", "R", colnames(normalized_counts))
)

dev.off()
log_message("Generated normalization effect plots")

# Save processed data
saveRDS(dds_filtered, "data/processed/dds_filtered_object.rds")
saveRDS(normalized_counts, "data/processed/normalized_counts.rds")

# Save filtering statistics
filtering_stats <- data.frame(
  total_genes = nrow(dds),
  genes_kept = sum(keep),
  genes_filtered = nrow(dds) - sum(keep),
  percentage_filtered = round((nrow(dds) - sum(keep))/nrow(dds) * 100, 2)
)

write.csv(filtering_stats, "results/tables/filtering_statistics.csv", row.names = FALSE)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "results/session_info_preprocessing.txt")

log_message("Preprocessing complete!")
