# PBMC Single-Cell RNA-seq Analysis

This repository provides a comprehensive analysis pipeline designed for the processing and visualization of single-cell RNA sequencing (scRNA-seq) data, specifically utilizing the Seurat package in the R programming environment. The pipeline is tailored to work with the PBMC 3K dataset, which is a widely used dataset provided by 10X Genomics. This dataset contains transcriptomic information derived from 3,000 individual Peripheral Blood Mononuclear Cells (PBMCs), offering a valuable resource for studying cellular heterogeneity within human blood. The repository includes various steps to process raw scRNA-seq data, from initial quality control to advanced visualizations, enabling users to gain insights into cellular compositions, gene expression patterns, and other key biological features at the single-cell level.

## Features
- Data Loading & Preprocessing: Reads 10X Genomics PBMC data and creates a Seurat object.
- Quality Control: Filters cells based on gene expression and mitochondrial gene content.
- Normalization & Scaling: Normalizes data, finds variable features, and scales expression values.
- Dimensionality Reduction: Uses PCA and UMAP for visualization.
- Clustering & Marker Gene Identification: Finds clusters and identifies top marker genes for each cluster.
- Visualization:
   - Violin plots (VlnPlot) for quality control and marker expression.
   - UMAP projection (DimPlot) to visualize clustering.
   - Feature expression plots (FeaturePlot) for specific marker genes.

## Plot Description
The included violin plots illustrate the expression distribution of key marker genes across identified cell clusters. Each violin plot shows the expression level of a gene (y-axis) across different clusters (x-axis). The width of each violin represents the density of expression values within each cluster.

## Output Files

- pbmc_tutorial.rds: Processed Seurat object for further analysis.
- Plots: Generated UMAP and violin plots for cell-type identification.

## Requirements

- Seurat
- ggplot2
- patchwork
- dplyr

## Usage
Run the script in an R environment to reproduce the analysis:

```r
source("pbmc_analysis.R")
```
