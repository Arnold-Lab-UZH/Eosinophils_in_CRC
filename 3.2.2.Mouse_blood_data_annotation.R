########## This code annotated cell type clusters from mouse blood data (non-tumor bearing mice) ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS(file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/tumor200_5000_25.rds")
#remove blood 
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("blood_wt"))

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(object = obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)
obj <- FindNeighbors(object = obj, dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, random.seed = 5, algorithm = 2)
obj <- RunUMAP(obj, dims = 1:15, seed.use = 5)
DimPlot(obj, label = TRUE)

##### Annotation
### automatic annotation using SingleR 
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj), ref = mouse.se, labels = mouse.se$label.main)
plotScoreHeatmap(results)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "seurat_clusters"
lapply(cell.types, function(x) project_annotation_to_umap(x, results, obj))

### DEGs
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "seurat_clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

##### rename clusters 
current.cluster.ids <- c(0:18)
new.cluster.ids <- c("Neutrophils","Neutrophils","Neutrophils","Neutrophils","Monocytes","Neutrophils","Neutrophils","Neutrophils",
                     "HSCs","HPCs","lowQ","Eosinophils","Monocytes","Monocytes","Basophils","T","Macrophages","NK","B")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj,reduction = "umap", label = TRUE, group.by = "annotation")

##### remove lowQ cells 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Neutrophils","Monocytes","HSCs","HPCs","Eosinophils","Basophils","T","Macrophages","NK","B"))

##### save object 
saveRDS(obj, "/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_blood_annotated.rds")
