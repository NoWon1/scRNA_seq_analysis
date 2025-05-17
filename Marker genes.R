library(Seurat)
library(dplyr)


mtx_file <- "Exp_data_UMIcounts.mtx"
features_file <- "Genes.tsv"
barcodes_file <- "barcodes.tsv"


sc_data <- ReadMtx(mtx = mtx_file, features = features_file,
                   cells = barcodes_file, feature.column = 1)


seurat_obj <- CreateSeuratObject(counts = sc_data)


seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")


VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                       nCount_RNA < 50000 & percent.mt < 10)


seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj)



seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)




seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)


DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "top10_cluster_markers.csv")


new.cluster.ids <- c(
  "CD4+ T cells",         
  "Skeletal muscle cells",
  "Macrophages",          
  "Epithelial cells",     
  "Secretory cells",      
  "Cytotoxic T/NK cells", 
  "Dendritic cells",      
  "Alveolar type II cells",
  "Fibroblasts",          
  "Adipocytes",           
  "B cells",              
  "Mast cells",           
  "Plasma cells",         
  "Cancer-associated cells",
  "Endothelial cells",    
  "Cycling cells (G2/M)", 
  "Cycling cells (S)",    
  "Cycling cells (G1)",   
  "Unknown1",             
  "Osteoclasts"           
)

names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)


DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.3)


DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()