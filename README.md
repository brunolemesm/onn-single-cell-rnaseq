# onn-single-cell-rnaseq
End-to-end single-cell RNA-seq analysis pipeline in R, integrating Seurat, Harmony, GPTCelltype (GPT-4-based cell type annotation), and CellChat for downstream communication and enrichment analyses.


# Single-cell RNA-seq Analysis Pipeline with Seurat, Harmony, and GPTCelltype

[![Made with R](https://img.shields.io/badge/Made%20with-R-blue?logo=r)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-%3E%3D5.0-orange)](https://satijalab.org/seurat/)
[![GPTCelltype](https://img.shields.io/badge/Cell%20Annotation-GPTCelltype-8A2BE2)](https://github.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This repository provides a comprehensive and reproducible workflow for single-cell RNA sequencing (scRNA-seq) analysis in **R**.  
The pipeline integrates **Seurat** and **Harmony** for data preprocessing, normalization, and clustering,  
uses **GPTCelltype** (a GPT-4–based annotation tool) for automated cell type prediction,  
and applies **CellChat** for intercellular communication inference.  

It includes:
- Quality control (QC) and filtering  
- Integration and batch correction with Harmony  
- Clustering, UMAP visualization, and marker detection  
- Automated annotation via GPTCelltype  
- Differential expression and enrichment analysis (GO and Reactome)  
- Cell–cell communication network analysis with CellChat  

**Goal:** Provide an automated, end-to-end workflow from raw 10X Genomics `.h5` files to biologically interpretable results.

## 🚀 Quick Start / How to Run

Follow these steps to run the single-cell RNA-seq pipeline locally.

### 1️⃣ Clone the repository

```bash
git clone https://github.com/seu-usuario/scRNAseq_pipeline_seurat_GPTCelltype.git
cd scRNAseq_pipeline_seurat_GPTCelltype
```

### 2️⃣ Install R and required packages

Make sure you have R >= 4.3 installed. Then install the required packages:
```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "tidyverse", "Matrix", "readxl", "openxlsx", "scales"))
BiocManager::install(c("SingleR", "org.Hs.eg.db", "ReactomePA", "clusterProfiler"))
install.packages(c("speckle", "EnhancedVolcano", "NMF", "ggalluvial", "pheatmap"))
# GPTCelltype should be installed according to its instructions
```

### 3️⃣ Prepare input data

Place your 10X Genomics filtered feature matrix .h5 files in a folder, e.g., data/.

Update the files vector in the R script to point to the correct paths:

```r

files <- c(
  "data/filtered_feature_bc_matrix_1.h5",
  "data/filtered_feature_bc_matrix_2.h5"
)
```
### 4️⃣ Run the pipeline

Open the main R script (scRNAseq_pipeline_seurat_GPTCelltype.R) and run it in R or RStudio:

```r

source("scRNAseq_pipeline_seurat_GPTCelltype.R")

```

### 5️⃣ Outputs

The pipeline generates:

QC plots (violin plots, barplots)

UMAPs by group, sample, and clusters

Differential expression results (xlsx)

Functional enrichment plots (GO, Reactome)

Cell-cell communication visualizations (CellChat)

DotPlots and heatmaps

All outputs are saved in the working directory or outputs/ folder (if created).
