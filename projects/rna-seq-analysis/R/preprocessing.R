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

# Create output directories if they don't exist
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/qc", recursive = TRUE, showWarnings = FALSE)

###################
# Data Download 
###################

# Download the dataset form GEO
gse <- getGEO("GSE286114", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))

# Get supplementary files
getGEOSuppFiles("GSE286114", baseDir = "data/raw")

###################
# Data Import 
###################

# Create temporary file and decompress
temp_file <- tempfile(fileext = ".xls")

if (file.exists("data/raw/GSE286114/GSE286114_gene_fpkm_RmpA_Escherichia.xls.gz")) {
  R.utils::gunzip("data/raw/GSE286114/GSE286114_gene_fpkm_RmpA_Escherichia.xls.gz", 
                  destname = temp_file, 
                  remove = FALSE)
} else {
  stop("File not found: Ensure GEO download was successful.")
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

# Save raw metadata
saveRDS(sample_metadata, "data/processed/sample_metadata.rds")

#Save gene information separately
gene_info <- counts_data %>% 
  select(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand, gene_length, gene_biotype, gene_description)

write_csv(gene_info, "data/processed/gene_info.csv")

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

# Create a histogram of the count distribution
hist(log10(total_counts_per_gene + 1), 
     breaks = 50, 
     main = "Distribution of Total Counts per Gene",
     xlab = "log10(Total Counts + 1)")


#Apply filtering 
#Filtering is applied to remove noise from our data which might affect the statistical analysis.

# Data-driven filtering
min_count <- 10  # minimum counts per sample
min_samples <- 3 # minimum number of samples
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds_filtered <- dds[keep,]


# Save filtering statistics
filtering_stats <- data.frame(
  total_genes = nrow(dds),
  genes_kept = sum(keep),
  genes_filtered = nrow(dds) - sum(keep)
)
print(filtering_stats)
write.csv(filtering_stats, "results/qc/filtering_statistics.csv")

# Visualize counts per sample before and after filtering
par(mfrow=c(1,2))
boxplot(log2(counts(dds) + 1), 
        main="Before Filtering",
        ylab="log2(counts + 1)",
        las=2)
boxplot(log2(counts(dds_filtered) + 1), 
        main="After Filtering",
        ylab="log2(counts + 1)",
        las=2)

###################
# QC Plots
###################

# Create QC plots for filtered data
pdf("figures/qc/filtered_data_qc.pdf")

# Sample correlation heatmap
vsd <- vst(dds_filtered, blind=TRUE)
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists),
         main="Sample Distance Matrix (Filtered Data)")

# Enhanced PCA plot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Filtered Data") +
  theme_minimal()

# Save filtered DESeq object
saveRDS(dds_filtered, "data/processed/dds_filtered_object.rds")

###################
# Normalization
###################
# Now proceed with normalization
dds_filtered <- estimateSizeFactors(dds_filtered)
normalized_counts <- counts(dds_filtered, normalized=TRUE)

# Save normalized data
saveRDS(normalized_counts, "data/processed/normalized_counts.rds")

# Create boxplots to visualize normalization effect
par(mfrow=c(1,2))
boxplot(log2(counts(dds_filtered)), 
        main="Raw counts", 
        ylab="log2(counts)", 
        las=2)
boxplot(log2(normalized_counts), 
        main="Normalized counts", 
        ylab="log2(counts)", 
        las=2)

saveRDS(normalized_counts, "data/processed/normalized_counts.rds")

# Save processing info
sessionInfo()
