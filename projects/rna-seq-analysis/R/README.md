# Analysis Pipeline Structure

## Files:
1. preprocessing.R
   - Raw data processing
   - Quality control
   - Data normalization
   - Saves processed data

2. analysis.Rmd
   - Loads preprocessed data
   - Creates visualizations
   - Performs differential expression analysis
   - Generates final report

## Workflow:
1. Run preprocessing.R first
2. Then knit analysis.Rmd

## Output Directory Structure:
- data/
  - raw/
  - processed/
- figures/
  - qc/
- results/
  - qc/
