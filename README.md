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

## Screenshots
![1](https://github.com/user-attachments/assets/47abdb4b-e26a-4985-8519-90d2fbba0f66)
![2](https://github.com/user-attachments/assets/a6b1dbab-3c15-4a21-8293-00fa9bc34b66)
![3](https://github.com/user-attachments/assets/5d837ba6-8ea8-44e6-9354-89bbf6370c96)
![4](https://github.com/user-attachments/assets/74b86de1-3a00-48bd-9189-0697fc72e42e)
![5](https://github.com/user-attachments/assets/77f6b3dc-8187-4809-a59b-307cf75dbad9)
![6](https://github.com/user-attachments/assets/82d9372b-c0a1-47a5-b97b-f906adfbd904)
