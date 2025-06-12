# Transcriptomics Workflows

A comprehensive comparison of R and Python approaches for RNA-seq differential expression analysis, demonstrating best practices in computational biology workflows.

## Overview

This repository showcases different computational approaches for analyzing the same RNA-seq dataset, comparing the strengths and applications of both R and Python ecosystems in transcriptomics research.

### Dataset: RmpA Regulation in E. coli
Analysis of transcriptional changes induced by RmpA (Regulator of Mucoid Phenotype A) overexpression in *Escherichia coli*, investigating its role in biofilm formation and metabolic regulation.

**Source**: GSE286114 (GEO Database)  
**Design**: 3 control vs 3 RmpA overexpression replicates  
**Platform**: RNA-sequencing  

## Workflows

### 🔵 [R Workflow](./r-workflow/)
**Status**: ✅ Complete

Comprehensive analysis using the Bioconductor ecosystem:
- **Preprocessing**: DESeq2, quality control
- **Differential Expression**: Statistical modeling with DESeq2
- **Pathway Analysis**: clusterProfiler, KEGG enrichment, GSEA
- **Visualization**: ggplot2, EnhancedVolcano, pheatmap

**Key Results**:
- 1,515 differentially expressed genes identified
- Ribosome biogenesis pathway strongly upregulated
- Metabolic reprogramming toward biofilm formation

[📊 View Complete Analysis](./r-workflow/analysis.html) | [📁 Explore R Code](./r-workflow/)

### 🐍 [Python Workflow](./python-workflow/)
**Status**: 🚧 In Development

Modern Python-based analysis pipeline:
- **Preprocessing**: pandas, scikit-learn
- **Differential Expression**: pyDESeq2, scipy.stats
- **Pathway Analysis**: gseapy, enrichr
- **Visualization**: matplotlib, seaborn, plotly

*Coming soon: Direct comparison of R and Python results*

### 📊 [Comparative Analysis](./comparative-analysis/)
**Status**: 📋 Planned

Head-to-head comparison of methodologies:
- Statistical concordance between R and Python approaches
- Performance benchmarking
- Workflow complexity analysis
- Reproducibility assessment

## Key Findings

| Aspect | Finding | Biological Significance |
|--------|---------|------------------------|
| **Ribosome Pathway** | 44 genes upregulated (p < 5e-11) | Enhanced protein synthesis capacity |
| **Nucleotide Sugars** | 15 genes upregulated (p < 3e-06) | Exopolysaccharide precursor production |
| **Fatty Acid Metabolism** | 10 genes downregulated (p < 4e-04) | Metabolic resource redirection |
| **Overall Impact** | 22% of transcriptome affected | Coordinated biofilm-associated reprogramming |

## Technologies Demonstrated

### R Ecosystem
![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-1f8c89?style=for-the-badge)

- DESeq2, clusterProfiler, EnhancedVolcano
- ggplot2, pheatmap, pathview
- R Markdown for reproducible reporting

### Python Ecosystem  
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Jupyter](https://img.shields.io/badge/Jupyter-F37626?style=for-the-badge&logo=jupyter&logoColor=white)

- pandas, numpy, scipy, scikit-learn
- matplotlib, seaborn, plotly
- Jupyter notebooks for interactive analysis

## Repository Structure

```
transcriptomics-workflows/
├── README.md                    # This overview
├── r-workflow/                  # Complete R-based analysis
│   ├── README.md
│   ├── analysis.Rmd
│   ├── preprocessing.R
│   ├── data/
│   ├── figures/
│   └── results/
├── python-workflow/             # Python implementation (in development)
│   ├── README.md
│   ├── analysis.ipynb
│   ├── preprocessing.py
│   └── requirements.txt
└── comparative-analysis/        # Cross-platform comparison
    └── method-comparison.md
```

## Quick Start

### R Workflow
```bash
git clone https://github.com/vikos77/transcriptomics-workflows.git
cd transcriptomics-workflows/r-workflow
Rscript preprocessing.R
# Open analysis.Rmd in RStudio
```

### Python Workflow (Coming Soon)
```bash
cd python-workflow
pip install -r requirements.txt
jupyter notebook analysis.ipynb
```

## Publications & References

**Dataset Source**: Yao S, Huang J, Geng J, et al. (2023). RmpA as a Global Regulator Modulates Switching Between Hypermucoviscosity and Biofilm. GEO: GSE286114.

**Methodological References**: 
- Love MI, et al. (2014). DESeq2. *Genome Biology*
- Yu G, et al. (2012). clusterProfiler. *OMICS*

## Contact

**Vigneshwaran Muthuraman**  
📧 vigneshwaran0594@gmail.com  
🔗 [LinkedIn](https://www.linkedin.com/in/vigneshwaran-muthuraman-24491746/)  
🐙 [GitHub](https://github.com/vikos77)

---

⭐ **Star this repository** if you find these workflows useful for your transcriptomics research!