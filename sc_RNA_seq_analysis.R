if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

pbmc.data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p1)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
print(p2)

print("Number of clusters:")
print(length(levels(Idents(pbmc))))

all.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
print("Top markers per cluster:")
print(top_markers)

features <- c("CCR7", "CD14", "IL7R", "CD79A", "CD8A", "FCGR3A", "NKG7", "FCER1A", "PPBP")
p3 <- FeaturePlot(pbmc, features = features, ncol = 3)
print(p3)

p4 <- VlnPlot(pbmc, features = features, ncol = 3)
print(p4)

cluster_num <- length(levels(Idents(pbmc)))
print(paste("Number of clusters to label:", cluster_num))

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

if(length(new.cluster.ids) != cluster_num) {
  stop(paste("Number of labels (", length(new.cluster.ids), ") doesn't match number of clusters (", cluster_num, ")", sep=""))
}

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "pbmc_tutorial.rds")
