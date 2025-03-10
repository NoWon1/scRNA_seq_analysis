#data link https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)
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

p1 <- VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "ctrl_qc_metrics.png"), p1, width = 12, height = 6)

p2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(paste0(out_path, "ctrl_feature_scatter.png"), p2, width = 12, height = 6)

ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
p3 <- VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "ctrl_filtered_metrics.png"), p3, width = 12, height = 6)

p4 <- VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "stim_qc_metrics.png"), p4, width = 12, height = 6)

p5 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(paste0(out_path, "stim_feature_scatter.png"), p5, width = 12, height = 6)

stim <- subset(stim, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)
p6 <- VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggsave(paste0(out_path, "stim_filtered_metrics.png"), p6, width = 12, height = 6)

ctrl <- NormalizeData(ctrl)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 3000)
stim <- NormalizeData(stim)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 3000)

immune.merged <- merge(ctrl, stim)
immune.merged <- NormalizeData(immune.merged)
immune.merged <- FindVariableFeatures(immune.merged, selection.method = "vst", nfeatures = 3000)
immune.merged <- ScaleData(immune.merged, vars.to.regress = c("percent.mt", "percent.rb"))
immune.merged <- RunPCA(immune.merged, npcs = 50)

p7 <- ElbowPlot(immune.merged, ndims = 50)
ggsave(paste0(out_path, "merged_elbow.png"), p7, width = 8, height = 6)

immune.integrated <- RunHarmony(immune.merged, "orig.ident")
p8 <- DimPlot(object = immune.integrated, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
ggsave(paste0(out_path, "harmony_dimplot.png"), p8, width = 8, height = 6)

immune.integrated <- FindNeighbors(immune.integrated, reduction = "harmony", dims = 1:30)
immune.integrated <- FindClusters(immune.integrated, resolution = c(0.4, 0.6, 0.8, 1.0))
immune.integrated <- RunUMAP(immune.integrated, reduction = "harmony", dims = 1:30)
immune.integrated <- RunTSNE(immune.integrated, reduction = "harmony", dims = 1:30)

p9 <- DimPlot(immune.integrated, reduction = "umap", group.by = "orig.ident", pt.size = .1)
p10 <- DimPlot(immune.integrated, reduction = "umap", group.by = "RNA_snn_res.1", label = TRUE, pt.size = .1)
p11 <- DimPlot(immune.integrated, reduction = "tsne", group.by = "orig.ident", pt.size = .1)
p12 <- DimPlot(immune.integrated, reduction = "tsne", group.by = "RNA_snn_res.1", label = TRUE, pt.size = .1)

combined_plot <- (p9 | p10) / (p11 | p12)
ggsave(paste0(out_path, "integrated_dimplots.png"), combined_plot, width = 14, height = 12)

DefaultAssay(immune.integrated) <- 'RNA'

marker_genes <- c("CD3D", "CD4", "CD8A", "GNLY", "NKG7", "FCGR3A", "MS4A1", "CD79A", "MZB1", "CST3", "LYZ", "IFI6", "ISG15", "MX1")
p13 <- FeaturePlot(immune.integrated, features = marker_genes, split.by = "orig.ident", 
                   max.cutoff = 3, cols = c("grey", "red"), reduction = "umap", ncol = 3)
ggsave(paste0(out_path, "marker_genes.png"), p13, width = 18, height = 20)

p14 <- DotPlot(immune.integrated, features = marker_genes, split.by = "orig.ident") + RotatedAxis()
ggsave(paste0(out_path, "marker_dotplot.png"), p14, width = 14, height = 8)

immune.cluster.markers <- FindAllMarkers(immune.integrated, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
write.csv(immune.cluster.markers, paste0(out_path, 'immune_cluster_markers.csv'))

top10 <- immune.cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p15 <- DoHeatmap(immune.integrated, features = top10$gene, label = TRUE)
ggsave(paste0(out_path, "top10_heatmap.png"), p15, width = 14, height = 20)

immune.integrated <- RenameIdents(immune.integrated, 
                                  `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                  `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", 
                                  `7` = "NK", `8` = "DC", `9` = "B Activated", `10` = "Mk", 
                                  `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

immune.integrated$CellTypes <- Idents(immune.integrated)
p16 <- DimPlot(immune.integrated, reduction = "umap", group.by = "CellTypes", pt.size = .1, label = TRUE, repel = TRUE)
ggsave(paste0(out_path, "final_celltypes.png"), p16, width = 10, height = 8)

ct_markers <- FindAllMarkers(immune.integrated, logfc.threshold = 0.25, only.pos = TRUE, min.pct = 0.25)
write.csv(ct_markers, paste0(out_path, 'celltype_markers.csv'))

cell_proportions <- table(Idents(immune.integrated), immune.integrated$orig.ident)
prop_table <- as.data.frame(prop.table(cell_proportions, margin = 2))
colnames(prop_table) <- c("CellType", "Condition", "Proportion")

p17 <- ggplot(prop_table, aes(x = Condition, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Type Proportions", x = "Condition", y = "Proportion")
ggsave(paste0(out_path, "celltype_proportions.png"), p17, width = 10, height = 8)

immune.integrated$celltype.stim <- paste(Idents(immune.integrated), immune.integrated$orig.ident, sep = "_")
Idents(immune.integrated) <- "celltype.stim"

DEG_list <- list()
cell_types <- levels(immune.integrated$CellTypes)
for(i in cell_types) {
  DEG_list[[i]] <- FindMarkers(immune.integrated, ident.1 = paste0(i, "_IMMUNE_STIM"), 
                               ident.2 = paste0(i, "_IMMUNE_CTRL"), 
                               logfc.threshold = 0.25, test.use = "wilcox")
  write.csv(DEG_list[[i]], paste0(out_path, 'DEG_', gsub(" ", "_", i), '.csv'))
}

aggregate_immune <- AggregateExpression(immune.integrated, 
                                        group.by = c("CellTypes", "orig.ident"), 
                                        return.seurat = TRUE)

interferon_genes <- c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p18 <- CellScatter(aggregate_immune, "CD14 Mono_IMMUNE_CTRL", "CD14 Mono_IMMUNE_STIM", highlight = interferon_genes)
p19 <- LabelPoints(plot = p18, points = interferon_genes, repel = TRUE)
ggsave(paste0(out_path, "monocyte_ifn_response.png"), p19, width = 10, height = 8)

saveRDS(immune.integrated, paste0(out_path, "immune_integrated_analysis.rds"))