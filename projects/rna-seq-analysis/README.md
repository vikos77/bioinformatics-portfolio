# Project: Transcriptional Analysis of RmpA-Mediated Regulation in Bacteria

## Introduction

This project investigates the transcriptional response of *Escherichia coli* to RmpA overexpression using RNA-Sequencing (RNA-Seq) data. RmpA (Regulator of Mucoid Phenotype A) is a bacterial protein known to play a role in [mention RmpA's known/suspected functions, e.g., virulence, capsule synthesis, biofilm formation in *Klebsiella pneumoniae*]. In this study, we analyze RNA-Seq data from *E. coli* cells with induced RmpA overexpression compared to control cells to identify genes and pathways regulated by RmpA in *E. coli*. This analysis aims to provide insights into the potential roles of RmpA in *E. coli*, particularly concerning its influence on bacterial physiology and biofilm-related processes.

## Key Findings

*   **Significant Differential Gene Expression:** Transcriptional profiling revealed significant changes in gene expression upon RmpA overexpression. We identified [Number of Up-regulated Genes - **FILL IN FROM ANALYSIS REPORT**] genes significantly up-regulated and [Number of Down-regulated Genes - **FILL IN FROM ANALYSIS REPORT**] genes significantly down-regulated (FDR < 0.05, |log2FoldChange| > 1).
*   **Ribosome Biogenesis Up-regulation:** Pathway analysis (KEGG ORA and GSEA) identified a strong enrichment of the Ribosome pathway, suggesting RmpA promotes increased protein synthesis.
*   **Nucleotide Sugar Biosynthesis Pathway Changes:**  The "Biosynthesis of various nucleotide sugars" pathway was also enriched, indicating a potential role for RmpA in regulating exopolysaccharide (EPS) production, a key component of bacterial biofilms.
*   **Metabolic Shifts:** Down-regulation of genes in 'Fatty acid degradation' and 'Glyoxylate and dicarboxylate metabolism' pathways suggests RmpA overexpression may induce metabolic reprogramming in *E. coli* cells.

**For a detailed analysis, including visualizations and complete pathway analysis results, please see the [Transcriptional Analysis Report](reports/analysis.html).**

## Preprocessing Analysis

The `R/preprocess_data.R` script performs the initial steps of the RNA-Seq data analysis, including:

*   **Data Download:**  Downloads the raw RNA-Seq dataset GSE286114 from the Gene Expression Omnibus (GEO).
*   **Data Import:** Imports the gene expression data (FPKM values) from the downloaded Excel file and creates a DESeq2 dataset object.
*   **Data Filtering:**
    *   Genes with very low counts across samples are filtered out to reduce noise and improve the power of differential expression analysis.
    *   **Filtering Rationale:** Low-count genes are often unreliable and can introduce noise into downstream analyses. Filtering helps to remove these noisy genes, focusing on genes with more robust and biologically meaningful expression changes.
    *   Filtering statistics (number of genes kept/filtered):

       ```
       total_genes genes_kept genes_filtered
     1        5013       4781            232
       ```

    *   Filtering statistics are saved to `results/qc/filtering_statistics.csv`.
    *   Boxplots visualizing count distributions **before and after filtering** are saved to `figures/qc/count_boxplots_filtering_effect.png` to **demonstrate the effect of filtering on data distribution.**

    **Visualization of Count Distribution After Filtering:**
    ![Boxplot of log2 counts after filtering](figures/qc/count_boxplots_filtering_effect_after.png)

   *   **Quality Control (QC):**
        *   **Sample Correlation Heatmap:**  A heatmap of sample-to-sample distances is generated to assess overall sample relationships and identify potential outliers. Saved to `figures/qc/filtered_data_qc.pdf`.
        *   **Principal Component Analysis (PCA):** PCA is performed to visualize sample clustering based on gene expression and check if samples group according to the experimental condition (control vs. RmpA overexpression). Saved to `figures/qc/filtered_data_qc.pdf`.

    **PCA Plot of Filtered Data:**
    ![PCA Plot of Filtered Data](figures/qc/filtered_data_qc.pdf)

   *   **Normalization:**
        *   Gene expression counts are normalized using the DESeq2 normalization method to account for differences in library size and sequencing depth between samples.
        *   Boxplots visualizing count distributions before and after normalization are saved to `figures/qc/count_boxplots_normalization_effect.pdf`.

    **Visualization of Count Distribution After Normalization:**

    ![Boxplots of Count Distribution After Normalization](figures/qc/count_boxplots_normalization_effect.png)

   *   **Data Saving:** Processed data objects (DESeq2 object, normalized counts, sample metadata, gene information) are saved as RDS files in the `data/processed/` directory for use in the main analysis report.
   *   **Session Information:**  The R session information (R version and package versions) is saved to `results/qc/session_info.txt` for reproducibility.


## Tools and Technologies Used

*   **R:** Statistical computing environment
*   **DESeq2:**  Differential gene expression analysis
*   **clusterProfiler:** KEGG pathway over-representation analysis and gene set enrichment analysis
*   **pathview:** KEGG pathway visualization
*   **pheatmap, EnhancedVolcano, ggplot2, tidyverse:** R packages for data visualization
*   **KEGGREST:**  KEGG database access in R
*   **org.EcK12.eg.db:** *E. coli* K12 annotation database

## Project Structure

rna-seq-analysis/
├── data/
│ └── processed/ # Processed RNA-Seq data (RDS objects, counts)
│ ├── dds_filtered_object.rds
│ └── normalized_counts.rds
│ └── raw/ # (Empty) Placeholder for raw data
│ └── README.md # Explanation of raw data storage
├── R/
│ └── analysis.Rmd # R Markdown report for RNA-Seq analysis
│ └── preprocess_data.R # R script for data preprocessing
├── reports/ # Generated reports
│ └── analysis.html # HTML report of RNA-Seq analysis
├── results/ # Analysis results summaries and outputs
│ └── qc/ # Quality control results
│ ├── filtering_statistics.csv
│ └── filtered_data_qc.pdf
│ └── count_boxplots_filtering_effect_before.png
│ └── count_boxplots_filtering_effect_after.png
│ └── count_boxplots_normalization_effect.pdf
│ └── pathway/ # Pathway analysis results
│ ├── kegg_gsea_results.csv
│ ├── kegg_ora_results.csv
│ └── pathway_enrichment_results.csv
└── README.md # Project overview and key findings