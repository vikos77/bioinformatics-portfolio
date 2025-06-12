#################################################
# RNA-seq Differential Expression Analysis
# Project: Fibrosis Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(ashr)
  library(org.Mm.eg.db)
})

# Logging function
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting differential expression analysis")

# Load DESeq2 object (assumes preprocessing has been done)
dds <- readRDS("results/dds_object.rds")

###################
# Run DESeq2 Analysis
###################

log_message("Running DESeq2 differential expression pipeline")

# Perform differential expression analysis with DESeq2
dds <- DESeq(dds)

# Verify result names
resultsNames(dds)

# Create output directories
dir.create("results/figures/differential_expression", recursive = TRUE, showWarnings = FALSE)

# Plot dispersion estimates
png("results/figures/differential_expression/dispersion_estimates.png", width=480, height=480)
plotDispEsts(dds, main="Dispersion Estimates")
dev.off()

###################
# Extract and Process Results
###################

log_message("Extracting differential expression results")

# Extract results for fibrosis vs. normal comparison
# alpha = 0.05 sets the significance cutoff for adjusted p-values
res <- results(dds,
               contrast = c("condition", "fibrosis", "normal"),
               alpha = 0.05
)

# Apply log fold change shrinkage to improve accuracy for visualization and ranking
log_message("Applying log fold change shrinkage")
res_shrunk <- lfcShrink(dds,
                        contrast = c("condition", "fibrosis", "normal"),
                        res = res,
                        type = "ashr")  # Using ashr for LFC shrinkage

# Display result summary
summary(res_shrunk)

# Create MA plots to visualize DE results before and after shrinkage
png("results/figures/differential_expression/ma_plots.png", width=480, height=480)
par(mfrow=c(1,2))
plotMA(res, main="Unshrunk LFC")
plotMA(res_shrunk, main="Shrunk LFC")
dev.off()

###################
# Annotate Results
###################

log_message("Annotating results with gene symbols")

# Convert results to a data frame for easier manipulation
res_df <- as.data.frame(res_shrunk) %>% 
  rownames_to_column("gene_id")

# Add gene symbols using the mouse genome annotation package
res_df$symbol <- mapIds(org.Mm.eg.db, 
                        keys = res_df$gene_id,
                        column = "SYMBOL", 
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Save the full results table
write.csv(res_df, "results/fibrosis_vs_normal_results.csv", row.names = FALSE)

###################
# Filter Significant Genes
###################

log_message("Filtering significant genes")

# Filter for significantly upregulated genes
up_genes <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange > 1)

# Filter for significantly downregulated genes
down_genes <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange < -1)

# Report the number of significant genes
log_message(paste("Significantly upregulated genes:", nrow(up_genes)))
log_message(paste("Significantly downregulated genes:", nrow(down_genes)))

# Save filtered gene lists
write.csv(up_genes, "results/upregulated_genes.csv", row.names = FALSE)
write.csv(down_genes, "results/downregulated_genes.csv", row.names = FALSE)

###################
# Visualization: Heatmap
###################

log_message("Creating expression heatmap for top DE genes")

# Load normalized counts
normalized_counts <- readRDS("results/normalized_counts.rds")

# Identify significant genes
sig_gene_ids <- res_df %>% filter(padj < 0.05) %>% pull(gene_id)
log_message(paste("Total significant genes:", length(sig_gene_ids)))

# If there are too many significant genes, take the top 50 by significance
if(length(sig_gene_ids) > 50) {
  sig_gene_ids <- res_df %>% 
    filter(padj < 0.05) %>% 
    arrange(padj) %>% 
    head(50) %>% 
    pull(gene_id)
}

# Extract normalized counts for significant genes
sig_norm_counts <- normalized_counts[sig_gene_ids,]

# Create heatmap
png("results/figures/differential_expression/expression_heatmap.png", width=480, height=480)
pheatmap(sig_norm_counts,
         color = colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         cluster_rows = TRUE,
         show_rownames = FALSE,
         annotation = data.frame(Condition = dds$condition,
                                 row.names = colnames(dds)),
         annotation_colors = list(Condition = c(fibrosis = "red", normal = "blue")),
         scale = "row",  # Z-score normalization by row
         main = "Top Differentially Expressed Genes")
dev.off()

###################
# Visualization: Volcano Plot
###################

log_message("Creating volcano plot")

# Add significance threshold indicators to result dataframe
res_df <- res_df %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1)

# Create enhanced volcano plot
png("results/figures/differential_expression/volcano_plot.png", width=480, height=480)
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red"), 
                     labels = c("Not significant", "Significant (padj < 0.05 & |log2FC| > 1)")) +
  labs(title = "Fibrosis vs Normal",
       subtitle = paste0(nrow(up_genes) + nrow(down_genes), " significantly differentially expressed genes"),
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
print(volcano_plot)
dev.off()

###################
# Visualization: Top Genes Expression Plot
###################

log_message("Creating top genes expression plot")

# Get top 10 differentially expressed genes by adjusted p-value
top10 <- res_df %>% 
  filter(padj < 0.05) %>% 
  arrange(padj) %>% 
  head(10) %>% 
  pull(gene_id)

# Extract normalized counts for top genes
top10_counts <- normalized_counts[top10, ]

# Create a tibble with gene symbol mapping for better plot labels
gene_labels <- res_df %>%
  filter(gene_id %in% top10) %>%
  select(gene_id, symbol) %>%
  mutate(symbol = ifelse(is.na(symbol), gene_id, symbol))  # Use Ensembl ID if no symbol

# Convert to long format for plotting
library(tidyr)
top10_long <- as.data.frame(top10_counts) %>%
  rownames_to_column("gene_id") %>%
  gather(key = "sample", value = "count", -gene_id) %>%
  inner_join(data.frame(sample = colnames(dds),
                        condition = dds$condition,
                        row.names = NULL),
             by = "sample") %>%
  # Add gene symbols for labeling
  left_join(gene_labels, by = "gene_id") %>%
  mutate(label = ifelse(is.na(symbol), gene_id, symbol))  # Use ID if no symbol

# Create enhanced expression plot
png("results/figures/differential_expression/top_genes_expression.png", width=480, height=480)
expression_plot <- ggplot(top10_long, aes(x = label, y = count, color = condition)) +
  geom_point(position = position_jitter(width = 0.2), size = 3, alpha = 0.8) +
  scale_y_log10() +  # Log scale helps visualize differences
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.title = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12)
  ) +
  labs(title = "Expression of Top 10 Differentially Expressed Genes",
       x = "Gene",
       y = "Normalized Count (log10 scale)",
       color = "Condition")
print(expression_plot)
dev.off()

log_message("Differential expression analysis complete")
