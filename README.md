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
