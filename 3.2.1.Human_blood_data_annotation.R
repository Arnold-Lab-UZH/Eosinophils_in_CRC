########## This code annotated cell type clusters from mouse blood data (non-tumor bearing mice) ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS(file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/blood_healthy_combined.rds")

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "condition",raster=FALSE)

obj <- JoinLayers(obj)

##### DEGs per cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =100, wt = avg_log2FC))

##### nFeature and percent.mito per cluster 
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")

##### SingleR for broad annotation
human.se <- celldex::NovershternHematopoieticData()
results <- SingleR(test = as.SingleCellExperiment(obj), ref = human.se, labels = human.se$label.main)
plotScoreHeatmap(results)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "mnn.clusters"
lapply(cell.types, function(x) project_annotation_to_umap_fastMNN(x, results, obj))

##### rename
current.cluster.ids <- c(0:12)
new.cluster.ids <- c("Eosinophils","Eosinophils","Eosinophils","Eosinophils",
                     "lowQ","Neutrophils","Eosinophils","other","Monocytes","lowQ","lowQ",
                     "T","other")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("Eosinophils","Neutrophils","T","Monocytes"))

##### DEGs per annotated cluster 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/eos_human/data_files/DEGs/DEGs_annotated_clusters_healthy_blood.csv")

##### save object 
saveRDS(obj, "/data/khandl/Eosinophils_in_CRC/seurat_objects/human_blood_annotated.rds")
