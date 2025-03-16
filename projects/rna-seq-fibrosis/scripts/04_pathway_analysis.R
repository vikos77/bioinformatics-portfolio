#################################################
# RNA-seq Pathway Enrichment Analysis
# Project: Fibrosis Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(enrichR)
  library(DOSE)
})

# Logging function
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting pathway enrichment analysis")

# Create output directories
dir.create("results/figures/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/enrichment", recursive = TRUE, showWarnings = FALSE)

# Load differential expression results
res_df <- read.csv("results/fibrosis_vs_normal_results.csv")
res <- readRDS("results/dds_object.rds")

###################
# Gene ID Conversion
###################

log_message("Converting gene IDs for enrichment analysis")

# Extract significant genes
sig_genes <- res_df %>% 
  filter(padj < 0.05) %>% 
  pull(gene_id)

# Convert Ensembl IDs to Entrez IDs (required for many enrichment tools)
entrez_ids <- bitr(sig_genes, 
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Mm.eg.db)

# Report mapping success rate
log_message(paste("Successfully mapped", nrow(entrez_ids), "out of", length(sig_genes), "genes"))

###################
# Separate Up and Down Regulated Genes
###################

log_message("Separating up and down regulated genes")

# Extract upregulated genes (padj < 0.05 & log2FC > 1)
up_genes <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange > 1) %>% 
  pull(gene_id)

# Extract downregulated genes (padj < 0.05 & log2FC < -1)
down_genes <- res_df %>% 
  filter(padj < 0.05 & log2FoldChange < -1) %>% 
  pull(gene_id)

# Convert to Entrez IDs
up_entrez <- bitr(up_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
  pull(ENTREZID)

down_entrez <- bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
  pull(ENTREZID)

# Report counts
log_message(paste("Upregulated genes:", length(up_genes), "(mapped:", length(up_entrez), ")"))
log_message(paste("Downregulated genes:", length(down_genes), "(mapped:", length(down_entrez), ")"))

###################
# GO Enrichment Analysis
###################

log_message("Performing GO enrichment analysis")

# Perform Gene Ontology enrichment analysis (Biological Process)
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "BP",  # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ego), "results/tables/enrichment/go_enrichment.csv", row.names = FALSE)

###################
# GO Visualization
###################

log_message("Creating GO enrichment visualizations")

# Barplot of enriched GO terms
png("results/figures/enrichment/go_barplot.png", width=480, height=480)
barplot(ego, showCategory=15) + 
  ggtitle("GO Biological Process Enrichment") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

# Dotplot of enriched GO terms
png("results/figures/enrichment/go_dotplot.png", width=480, height=480)
dotplot(ego, showCategory=15) + 
  ggtitle("GO Enrichment Dot Plot") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

# Enrichment map - network of related terms
png("results/figures/enrichment/go_network.png", width=480, height=480)
emapplot(pairwise_termsim(ego), showCategory = 20) +
  ggtitle("Network of Enriched GO Terms") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

###################
# KEGG Pathway Analysis
###################

log_message("Performing KEGG pathway analysis")

# Perform KEGG pathway enrichment analysis
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                    organism = 'mmu',  # Mouse
                    pvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ekegg), "results/tables/enrichment/kegg_enrichment.csv", row.names = FALSE)

# Barplot of enriched KEGG pathways
png("results/figures/enrichment/kegg_barplot.png", width=480, height=480)
barplot(ekegg, showCategory=15) + 
  ggtitle("KEGG Pathway Enrichment") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

###################
# Compare Up vs Down Regulated Pathways
###################

log_message("Comparing enrichment in up vs down regulated genes")

# Create a named list for compareCluster
gene_list <- list(Up = up_entrez, Down = down_entrez)

# Compare GO terms between up and down regulated genes
compare_go <- compareCluster(gene_list, 
                             fun = "enrichGO",
                             OrgDb = org.Mm.eg.db,
                             ont = "BP")

# Save comparison results
write.csv(as.data.frame(compare_go), "results/tables/enrichment/up_vs_down_go.csv", row.names = FALSE)

# Visualize GO comparison
png("results/figures/enrichment/up_vs_down_go_dotplot.png", width=480, height=480)
dotplot(compare_go, showCategory = 10) + 
  ggtitle("Up vs Down Regulated Genes: GO Enrichment") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

# Compare KEGG pathways
compare_kegg <- compareCluster(gene_list,
                               fun = "enrichKEGG",
                               organism = "mmu")

# Save KEGG comparison results
write.csv(as.data.frame(compare_kegg), "results/tables/enrichment/up_vs_down_kegg.csv", row.names = FALSE)

# Visualize KEGG comparison
png("results/figures/enrichment/up_vs_down_kegg_dotplot.png", width=480, height=480)
dotplot(compare_kegg, showCategory = 10) + 
  ggtitle("Up vs Down Regulated Genes: KEGG Pathways") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

###################
# EnrichR Analysis
###################

log_message("Performing EnrichR analysis for additional databases")

# Get gene symbols for EnrichR
gene_symbols <- bitr(sig_genes, 
                     fromType = "ENSEMBL", 
                     toType = "SYMBOL", 
                     OrgDb = org.Mm.eg.db) %>%
  pull(SYMBOL)

# Select databases relevant for mouse studies
selected_dbs <- c("GO_Biological_Process_2021", 
                  "GO_Molecular_Function_2021",
                  "GO_Cellular_Component_2021",
                  "KEGG_2021_Human",  # Many mouse genes have human orthologs
                  "WikiPathways_2019_Mouse",
                  "Mouse_Gene_Atlas")

# Run enrichment
enriched <- enrichr(gene_symbols, selected_dbs)

# Save EnrichR results
for(db in names(enriched)) {
  write.csv(enriched[[db]], 
            file = paste0("results/tables/enrichment/enrichr_", gsub("[^a-zA-Z0-9]", "_", db), ".csv"), 
            row.names = FALSE)
}

# Plot EnrichR results for GO Biological Process
png("results/figures/enrichment/enrichr_plot.png", width=480, height=480)
plotEnrich(enriched[["GO_Biological_Process_2021"]], 
           showTerms = 20, 
           numChar = 50, 
           y = "Count", 
           orderBy = "P.value") +
  ggtitle("Enrichment analysis by EnrichR") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
dev.off()

###################
# Gene Set Enrichment Analysis (GSEA)
###################

log_message("Performing Gene Set Enrichment Analysis (GSEA)")

# Prepare ranked gene list for GSEA
log_message("Preparing ranked gene list for GSEA")

# Extract results for fibrosis vs. normal comparison
res <- results(dds,
               contrast = c("condition", "fibrosis", "normal"),
               alpha =0.05
)
# Get all results with p-values
all_results <- as.data.frame(res)

# Keep only genes that can be mapped to Entrez IDs
mapped_genes <- bitr(rownames(all_results), 
                     fromType = "ENSEMBL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)

# Merge results with mappings
gene_data <- all_results %>%
  rownames_to_column("ENSEMBL") %>%
  inner_join(mapped_genes, by = "ENSEMBL")

# Handle duplicates - keep entry with largest absolute fold change
gene_data <- gene_data %>%
  group_by(ENTREZID) %>%
  slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

# Create ranking for GSEA (log2FoldChange)
gene_list <- gene_data$log2FoldChange
names(gene_list) <- gene_data$ENTREZID

# Sort the list for GSEA
gene_list <- sort(gene_list, decreasing = TRUE)

# Final checks and cleanup
gene_list <- na.omit(gene_list)
gene_list <- gene_list[is.finite(gene_list)]
if(any(duplicated(names(gene_list)))) {
  warning("Duplicate gene IDs found. Keeping first instance.")
  gene_list <- gene_list[!duplicated(names(gene_list))]
}

# Run GSEA with the cleaned data
gsea_result <- gseGO(geneList = gene_list,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# Save GSEA results
write.csv(as.data.frame(gsea_result), "results/tables/enrichment/gsea_results.csv", row.names = FALSE)

# Visualize GSEA results
png("results/figures/enrichment/gsea_plot.png", width=480, height=480)
gseaplot2(gsea_result, geneSetID = 1:3, title = "GSEA for Top Pathways")
dev.off()

# Create additional GSEA visualizations
if(length(gsea_result$ID) > 0) {
  # Ridge plot for multiple pathways
  png("results/figures/enrichment/gsea_ridgeplot.png", width=480, height=480)
  ridgeplot(gsea_result, showCategory = 10) + 
    ggtitle("GSEA Pathway Ridge Plot") +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  dev.off()
}

log_message("Pathway enrichment analysis complete")
