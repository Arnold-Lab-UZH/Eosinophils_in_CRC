########## This code integrated human eosinophils from tumor, NAT and blood healthy ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
blood <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_blood_annotated.rds")
colon <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")

# extract eosinophils 
Idents(blood) <- "annotation"
blood_eos <- subset(blood, idents = "Eosinophils")
Idents(colon) <- "annotation"
colon_eos <- subset(colon, idents = "Eosinophils")

obj <-  merge(blood_eos, y = c(colon_eos),
            add.cell.ids = c("blood", "colon"))
obj <- JoinLayers(obj)

##### pre-processing  
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$tissue)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE)

## check quality of clusters
VlnPlot(obj, features = "nFeature_RNA")
VlnPlot(obj, features = "percent.mt")
#remove cluster 9 (very high mito and very low nFeatures)
Idents(obj) <- "mnn.clusters"
obj <- subset(obj, idents = c(0,1,2,3,4,5,6,7,8))
obj <- JoinLayers(obj)

##### DEGs per cluster 
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

#cluster 8 has very high features and expressed also macrophages (CCR2) and mast cell/basophils (MS4A2,CPA3), DCs (FCER1A), neutropils (MS4A3), genes 
#--> remove it as doublets/potential mast/basophils
obj <- subset(obj, idents = c(0,1,2,3,4,5,6,7))
obj <- JoinLayers(obj)
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE)

##### combine clusters that are transcriptionally very similar 
current.cluster.ids <- c(0:7)
new.cluster.ids <- c(1,2,1,1,1,1,1,1)
obj$mnn.clusters <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label=TRUE) 

p <- DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE, label = TRUE, cols = c("#CA13EA","#EAA10F"))
ggsave(file = "/scratch/khandl/eos_human/eos_clustering/umap_clusters.svg", width = 8, height = 6, plot = p)

##### eosinophil markers 
markers <- c("CLC","CCR3","SYNE1","ADGRE5","SIGLEC8","IL5RA")
markers <- c("MARCHF3","FFAR2","MCTP2","RELB","CCL4L2","DACH1","CCL4L2","SSH2","IL1RL1","SYNE1","CEACAM1")
for(i in markers) {
  p <- FeaturePlot(obj, features = i, reduction = "umap.mnn", pt.size = 0.1) + scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,5))
  ggsave(paste0("/scratch/khandl/eos_human/eos_clustering/",i,".pdf"), width = 8, height = 5, plot = p)
} 

##### DEGs per cluster 
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/eos_human/eos_clustering/DEGs_per_cluster.csv")

## plot gene of interest per cluster 
markers_to_plot <- c("SMPD3","IL5RA","PRSS33","PIK3R6","ITGAL","ITGA4","ALOX5","ALOX15","TGFBR2","TGFB1", "RARA", #cluster 1 specific 
                     "CCL3","CCL3L1","CCL4","CCL4L2","NFKB1","NFKB2","REL", "RELB",
                     "TNFAIP2","TNFAIP3","CD69","AREG","HLA-DRA" ,"FOSB" #cluster 2 specific
                     )

p <- DotPlot(obj, features = markers_to_plot , dot.scale = 10, scale = FALSE) + RotatedAxis()
ggsave(file = "/scratch/khandl/eos_human/eos_clustering/dotPlot_genes_of_interest.svg", width = 8, height = 6, plot = p)

##### proportion of mnn.clusters per tissue
create_table_cell_type_prop(obj, "tissue","mnn.clusters","/scratch/khandl/eos_human/eos_clustering/","clusters")
df <- read.csv("/scratch/khandl/eos_human/eos_clustering/clusters_proportions_tissue_mnn.clusters.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#CA13EA","#EAA10F")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/eos_clustering/prop_barplot_clusters.svg", width = 12, height = 6, plot = p)

##### statistics between NAT and tumor
#first control, then the condition
Idents(obj) <- "mnn.clusters"
obj <- subset(obj, idents = c("1","2"))
cell_type_prop_stats(obj,"mnn.clusters","tissue_ctrl","tumor","tissue",1.41,
                     "/scratch/khandl/eos_human/eos_clustering/stats_tissue_ctrl_vs_tumor.svg") 

##### QC 
### Vlnplot of quality measures per cluster 
Idents(obj) <- "mnn.clusters"
p <- VlnPlot(obj, features = c("nFeature_RNA"), pt.size = 0, cols = c("#CA13EA","#EAA10F"))
ggsave("/scratch/khandl/eos_human/QC_cluster/nFeature.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(obj, features = c("percent.mt"), pt.size = 0 ,cols = c("#CA13EA","#EAA10F"))
ggsave("/scratch/khandl/eos_human/QC_cluster/mito.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(obj, features = c("nCount_RNA"), pt.size = 0, cols = c("#CA13EA","#EAA10F"))
ggsave("/scratch/khandl/eos_human/QC_cluster/nCount.svg", width = 8, height = 5, plot = p)

Idents(obj) <- "tissue"
p <- VlnPlot(obj, features = c("nFeature_RNA"), pt.size = 0,cols = c("#4168B2","#E91F44","#F4E800"))
ggsave("/scratch/khandl/eos_human/QC_tissues/nFeature.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(obj, features = c("percent.mt"), pt.size = 0,cols = c("#4168B2","#E91F44","#F4E800"))
ggsave("/scratch/khandl/eos_human/QC_tissues/mito.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(obj, features = c("nCount_RNA"), pt.size = 0,cols = c("#4168B2","#E91F44","#F4E800"))
ggsave("/scratch/khandl/eos_human/QC_tissues/nCount.svg", width = 8, height = 5, plot = p)

### umaps per sample 
#per tissue 
p <- DimPlot(obj,reduction = "umap.mnn",group.by = "tissue",raster=FALSE, cols = c( "#0667F4","#F4E806","#F4063F"))
ggsave(file = "/scratch/khandl/eos_human/eos_clustering/umap_tisues.svg", width = 8, height = 6, plot = p)

# blood samples 
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "blood_healhty")
p <- DimPlot(sub,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
ggsave(file = "/scratch/khandl/eos_human/eos_clustering/umap_blood_per_sample.svg", width = 8, height = 6, plot = p)

# tumor
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "tumor")
p <- DimPlot(sub,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
ggsave(file = "/scratch/khandl/eos_human/eos_clustering/umap_tumor_per_sample.svg", width = 8, height = 6, plot = p)

# NAT
Idents(obj) <- "tissue"
sub <- subset(obj, idents = "tissue_ctrl")
p <- DimPlot(sub,reduction = "umap.mnn",group.by = "condition",raster=FALSE)
ggsave(file = "/scratch/khandl/eos_human/eos_clustering/umap_tissue_per_sample.svg", width = 8, height = 6, plot = p)

##### save object 
saveRDS(obj, "/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated.rds")

