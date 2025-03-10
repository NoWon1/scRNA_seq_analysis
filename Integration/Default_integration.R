#data link https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data
ctrl.data <- read.table(file = "C:/Users/adity/OneDrive/Desktop/Coding/R/Integration/immune_alignment_expression_matrices/immune_control_expression_matrix.txt.gz", sep = "\t", header = TRUE, row.names = 1)
stim.data <- read.table(file = "C:/Users/adity/OneDrive/Desktop/Coding/R/Integration/immune_alignment_expression_matrices/immune_stimulated_expression_matrix.txt.gz", sep = "\t", header = TRUE, row.names = 1)

# Create Seurat objects for control data
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 3)
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")

# Visualize QC metrics for control
p1 <- VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p1)

p2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
      geom_smooth(method = "lm")
print(p2)

# Filter control cells
ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
p3 <- VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p3)

# Create Seurat objects for stimulated data
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 3)
stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")

# Visualize QC metrics for stimulated
p4 <- VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p4)

p5 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
      geom_smooth(method = "lm")
print(p5)

# Filter stimulated cells
stim <- subset(stim, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
p6 <- VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p6)

# Merge datasets for initial analysis
immune.merged <- merge(ctrl, stim)

# Process merged data
immune.merged <- NormalizeData(immune.merged, normalization.method = "LogNormalize", scale.factor = 10000)
immune.merged <- FindVariableFeatures(immune.merged, selection.method = "vst", nfeatures = 2000)

# Visualize variable features
top10 <- head(VariableFeatures(immune.merged), 10)
plot1 <- VariableFeaturePlot(immune.merged)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1 + plot2)

# Scale data and run PCA
all.genes <- rownames(immune.merged)
immune.merged <- ScaleData(immune.merged, features = all.genes)
immune.merged <- RunPCA(immune.merged, npcs = 30, verbose = FALSE)

# Determine dimensionality
ElbowPlot(immune.merged, ndims = 30)

# Cluster and visualize merged data
immune.merged <- FindNeighbors(immune.merged, reduction = "pca", dims = 1:30)
immune.merged <- RunUMAP(immune.merged, reduction = "pca", dims = 1:30)
immune.merged <- FindClusters(immune.merged, resolution = 0.5)

# Visualization of merged data
p7 <- DimPlot(immune.merged)
print(p7)

p8 <- DimPlot(immune.merged, group.by = 'orig.ident')
print(p8)

p9 <- DimPlot(immune.merged, split.by = 'orig.ident')
print(p9)

p10 <- DimPlot(immune.merged, split.by = 'orig.ident', label = TRUE) + NoLegend()
print(p10)

p11 <- FeaturePlot(immune.merged, features = "GNLY")
print(p11)

# Properly prepare individual datasets for integration
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Select integration features with fewer features
immune.integrated.features <- SelectIntegrationFeatures(object.list = list(ctrl, stim), nfeatures = 2000)

# Set memory limit
options(future.globals.maxSize = 200000000000)

# Run integration with proper parameters
immune.integrated.anchors <- FindIntegrationAnchors(
  object.list = list(ctrl, stim),
  anchor.features = immune.integrated.features, 
  reduction = "cca",
  dims = 1:20,
  seed = TRUE,
  verbose = FALSE
)

# Continue with integration
immune.integrated <- IntegrateData(anchorset = immune.integrated.anchors, verbose = FALSE)

# Process integrated data
DefaultAssay(immune.integrated) <- "integrated"
immune.integrated <- ScaleData(immune.integrated)
immune.integrated <- RunPCA(immune.integrated, npcs = 30, verbose = FALSE)
immune.integrated <- RunUMAP(immune.integrated, dims = 1:30)
immune.integrated <- FindNeighbors(immune.integrated, dims = 1:30)
immune.integrated <- FindClusters(immune.integrated, resolution = 0.5)

# Visualize integrated data
p12 <- DimPlot(immune.integrated, group.by = 'orig.ident')
print(p12)

p13 <- DimPlot(immune.integrated, split.by = 'orig.ident')
print(p13)

p14 <- DimPlot(immune.integrated, split.by = 'orig.ident', label = TRUE) + NoLegend()
print(p14)

# Switch to RNA assay for feature plots
DefaultAssay(immune.integrated) <- 'RNA'
p15 <- FeaturePlot(immune.integrated, features = c("CD3D", "GNLY", "IFI6"), 
            split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red"), reduction = "umap")
print(p15)

p16 <- FeaturePlot(immune.integrated, features = "GNLY")
print(p16)

# Find and save marker genes
out_path <- "D:/sc_analysis/"
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

immune.cluster.markers <- FindAllMarkers(immune.integrated, logfc.threshold = 0.1, only.pos = TRUE)
write.csv(immune.cluster.markers, paste(out_path, 'immune_cluster_markers.csv', sep = ''))

# Visualize top markers
top5_markers <- immune.cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(immune.integrated, features = top5_markers$gene, size = 3, angle = 90) + NoLegend()

# Rename cell types
immune.integrated <- RenameIdents(immune.integrated, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                  `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
                                  `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

# Store cell type identities in metadata
immune.integrated$cell_type <- Idents(immune.integrated)
p17 <- DimPlot(immune.integrated, label = TRUE, repel = TRUE) + NoLegend()
print(p17)

# Find cell type markers
ct_markers <- FindAllMarkers(immune.integrated, logfc.threshold = 0.1, only.pos = TRUE)

# Violin plots for key markers
VlnPlot(immune.integrated, features = c("CD3D", "CD14", "FCER1A", "FCGR3A", "CD8A", "MS4A1"), 
        pt.size = 0, ncol = 3)

# Dot plot for markers
DotPlot(immune.integrated, features = c("CD3D", "CD14", "FCER1A", "FCGR3A", "CD8A", "MS4A1"), 
        dot.scale = 8) + RotatedAxis()

# t-SNE visualization
immune.integrated <- RunTSNE(immune.integrated, dims = 1:30)
DimPlot(immune.integrated, reduction = "tsne", label = TRUE, repel = TRUE)

# Save integrated object
saveRDS(immune.integrated, file = paste0(out_path, "immune_integrated.rds"))

# Cell type proportion analysis
cell_counts <- table(Idents(immune.integrated), immune.integrated$orig.ident)
cell_props <- prop.table(cell_counts, margin = 2)
cell_props_df <- as.data.frame(cell_props)
colnames(cell_props_df) <- c("CellType", "Condition", "Proportion")

ggplot(cell_props_df, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Condition")

# Differential expression analysis
differential_expression <- FindMarkers(immune.integrated, ident.1 = "IMMUNE_STIM", 
                                      ident.2 = "IMMUNE_CTRL", 
                                      group.by = "orig.ident",
                                      subset.ident = "CD14 Mono")

# Install EnhancedVolcano if needed
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("EnhancedVolcano")
}

# Volcano plot
library(EnhancedVolcano)
EnhancedVolcano(differential_expression, 
                lab = rownames(differential_expression),
                x = "avg_log2FC",
                y = "p_val_adj",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0)
