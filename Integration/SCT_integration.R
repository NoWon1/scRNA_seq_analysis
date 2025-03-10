#data link https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)
library(cowplot)

plan("multiprocess", workers = 4)

set.seed(42)
out_path <- "D:/sc_analysis/"
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

ctrl.data <- read.table(file = "immune_alignment_expression_matrices/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "immune_alignment_expression_matrices/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 3, min.features = 200)
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 3, min.features = 200)

ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")

ctrl[["percent.rb"]] <- PercentageFeatureSet(ctrl, pattern = "^RP[SL]")
stim[["percent.rb"]] <- PercentageFeatureSet(stim, pattern = "^RP[SL]")

ctrl.metrics <- VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "ctrl_qc_metrics.png"), ctrl.metrics, width = 12, height = 6)

stim.metrics <- VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "stim_qc_metrics.png"), stim.metrics, width = 12, height = 6)

plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("Control")
plot2 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("Stimulated")
combined_scatter <- plot1 + plot2
ggsave(paste0(out_path, "feature_scatter.png"), combined_scatter, width = 12, height = 6)

ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
stim <- subset(stim, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)

ctrl.filtered <- VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
stim.filtered <- VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "filtered_metrics.png"), ctrl.filtered / stim.filtered, width = 12, height = 10)

ctrl <- NormalizeData(ctrl)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
ctrl <- ScaleData(ctrl, features = rownames(ctrl))
ctrl <- RunPCA(ctrl, features = VariableFeatures(ctrl))
ElbowPlot(ctrl, ndims = 50)
ggsave(paste0(out_path, "ctrl_elbow.png"), width = 8, height = 6)

ctrl <- SCTransform(ctrl, vars.to.regress = c("percent.mt", "percent.rb"))

stim <- NormalizeData(stim)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)
stim <- ScaleData(stim, features = rownames(stim))
stim <- RunPCA(stim, features = VariableFeatures(stim))
ElbowPlot(stim, ndims = 50)
ggsave(paste0(out_path, "stim_elbow.png"), width = 8, height = 6)

stim <- SCTransform(stim, vars.to.regress = c("percent.mt", "percent.rb"))

immune.features <- SelectIntegrationFeatures(object.list = c(ctrl, stim), nfeatures = 3000)
options(future.globals.maxSize = 200000000000)
immune.list <- PrepSCTIntegration(object.list = c(ctrl, stim), anchor.features = immune.features, verbose = FALSE)

immune.anchors <- FindIntegrationAnchors(object.list = immune.list, normalization.method = "SCT", 
                                         anchor.features = immune.features, verbose = FALSE)

immune.integrated <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", verbose = FALSE)

immune.integrated <- RunPCA(immune.integrated, verbose = FALSE, npcs = 50)
immune.integrated <- FindNeighbors(immune.integrated, dims = 1:50)
immune.integrated <- FindClusters(immune.integrated, resolution = c(0.3, 0.5, 0.8, 1.0))
immune.integrated <- RunUMAP(immune.integrated, dims = 1:50)
immune.integrated <- RunTSNE(immune.integrated, dims = 1:50)

p1 <- DimPlot(immune.integrated, reduction = "umap", group.by = 'orig.ident')
p2 <- DimPlot(immune.integrated, reduction = "tsne", group.by = 'orig.ident')
combined_plot <- p1 + p2
ggsave(paste0(out_path, "integrated_dimplot.png"), combined_plot, width = 12, height = 6)

p3 <- DimPlot(immune.integrated, reduction = "umap", split.by = 'orig.ident')
p4 <- DimPlot(immune.integrated, reduction = "tsne", split.by = 'orig.ident')
ggsave(paste0(out_path, "split_dimplot.png"), p3, width = 12, height = 6)
ggsave(paste0(out_path, "split_tsne.png"), p4, width = 12, height = 6)

labeled_clusters <- DimPlot(immune.integrated, reduction = "umap", split.by = 'orig.ident', label = TRUE) + NoLegend()
ggsave(paste0(out_path, "labeled_clusters.png"), labeled_clusters, width = 12, height = 6)

DefaultAssay(immune.integrated) <- 'SCT'

marker_genes <- c("CD3D", "CD4", "CD8A", "GNLY", "NKG7", "FCGR3A", "MS4A1", "CD79A", "MZB1", "CST3", "LYZ", "PPBP", "IFI6")
feature_plot <- FeaturePlot(immune.integrated, features = marker_genes, 
                            split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red"), reduction = "umap", ncol = 3)
ggsave(paste0(out_path, "marker_genes.png"), feature_plot, width = 18, height = 20)

dot_plot <- DotPlot(immune.integrated, features = marker_genes, split.by = "orig.ident") + RotatedAxis()
ggsave(paste0(out_path, "marker_dotplot.png"), dot_plot, width = 14, height = 8)

immune.cluster.markers <- FindAllMarkers(immune.integrated, logfc.threshold = 0.1, only.pos = TRUE, min.pct = 0.25)
write.csv(immune.cluster.markers, paste0(out_path, 'immune_cluster_markers.csv'))

top10 <- immune.cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heat_map <- DoHeatmap(immune.integrated, features = top10$gene, label = TRUE) + NoLegend()
ggsave(paste0(out_path, "top10_heatmap.png"), heat_map, width = 14, height = 20)

immune.integrated <- RenameIdents(immune.integrated, 
                                  `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                  `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", 
                                  `7` = "NK", `8` = "DC", `9` = "B Activated", `10` = "Mk", 
                                  `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

immune.integrated[["cell_type"]] <- Idents(immune.integrated)
final_plot <- DimPlot(immune.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
ggsave(paste0(out_path, "final_celltypes.png"), final_plot, width = 10, height = 8)

ct_markers <- FindAllMarkers(immune.integrated, logfc.threshold = 0.1, only.pos = TRUE, min.pct = 0.25)
write.csv(ct_markers, paste0(out_path, 'celltype_markers.csv'))

saveRDS(immune.integrated, paste0(out_path, "immune_integrated_analysis.rds"))

cell_proportions <- table(Idents(immune.integrated), immune.integrated$orig.ident)
prop_table <- as.data.frame(prop.table(cell_proportions, margin = 2))
colnames(prop_table) <- c("CellType", "Condition", "Proportion")

prop_plot <- ggplot(prop_table, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Type Proportions", x = "Condition", y = "Proportion")
ggsave(paste0(out_path, "celltype_proportions.png"), prop_plot, width = 10, height = 8)

DEG_list <- list()
cell_types <- levels(Idents(immune.integrated))
for(i in cell_types) {
  DEG_list[[i]] <- FindMarkers(immune.integrated, ident.1 = "IMMUNE_STIM", ident.2 = "IMMUNE_CTRL", 
                               group.by = "orig.ident", subset.ident = i, 
                               logfc.threshold = 0.25, test.use = "wilcox")
  write.csv(DEG_list[[i]], paste0(out_path, 'DEG_', gsub(" ", "_", i), '.csv'))
}