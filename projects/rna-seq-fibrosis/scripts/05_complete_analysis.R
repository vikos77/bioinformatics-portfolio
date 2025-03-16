#################################################
# Complete RNA-seq Analysis Pipeline
# Project: Fibrosis Analysis
# Author: Vigneshwaran Muthuraman
#################################################

# This script runs all analysis steps in sequence
# It's useful for reproducing the entire analysis workflow

# Initialize logging
cat("Starting complete RNA-seq analysis workflow\n")

# Create all necessary directories
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Run each script in sequence
cat("Step 1: Data preprocessing\n")
source("scripts/01_data_preprocessing.R")

cat("Step 2: Exploratory analysis\n")
source("scripts/02_exploratory_analysis.R")

cat("Step 3: Differential expression analysis\n")
source("scripts/03_differential_expression.R")

cat("Step 4: Pathway analysis\n")
source("scripts/04_pathway_analysis.R")

cat("Analysis complete. Results are available in the 'results' directory.\n")