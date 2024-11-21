########## This code integrates colon and tumor data from mouse samples and annotates clusters ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS(file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/tumor200_5000_25.rds")
#remove blood 
Idents(obj) <- "condition"
obj <- subset(obj, idents = c("tumor_wt","disseminated_wt","adult_colon_wt", "tumor_phil", "disseminated_phil", "adult_colon_phil"))

##### pre-processing
obj <- NormalizeData(obj,normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(object = obj, features = VariableFeatures(object =obj), npcs = 20, verbose = FALSE)

##### fastMNN integration 
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$condition)
obj <- IntegrateLayers(object = obj, method = FastMNNIntegration,new.reduction = "integrated.mnn",
                       verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:15)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn",return.model = TRUE)
DimPlot(obj,reduction = "umap.mnn", label = TRUE)

obj <- JoinLayers(obj)

##### Broad cell type annotation using SingleR 
mouse.se <- celldex::ImmGenData()
results <- SingleR(test = as.SingleCellExperiment(obj), ref = mouse.se, labels = mouse.se$label.main)
plotScoreHeatmap(results)
cell.types <- unique(results$pruned.labels)
Idents(obj) <- "mnn.clusters"
lapply(cell.types, function(x) project_annotation_to_umap_fastMNN(x, results, obj))

### DEGs
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
DimPlot(obj,reduction = "umap.mnn",raster=FALSE, group.by = "mnn.clusters",label = TRUE)
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

### rename clusters 
current.cluster.ids <- c(0:18)
new.cluster.ids <- c("Neutrophils","B","Mono_Mac","ILC_T","B","Mono_Mac","Mono_Mac","Eosinophils","DCs","lowQ","Mono_Mac",
                     "ILC_T","Neutrophils","B","B","ILC_T","Neutrophils","ILC_T","B")
obj$annotation <- plyr::mapvalues(x = obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj,reduction = "umap.mnn", label = TRUE, group.by = "annotation1")

##### Sub-clustering 
### Neutrophils 
Idents(obj) <- "annotation"
subCl <- FindSubCluster(obj,cluster = "Neutrophils",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.5)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Neutrophils_0","Neutrophils_1","Neutrophils_2","Neutrophils_3","Neutrophils_4",
                                         "Neutrophils_5","Neutrophils_6"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

## rename
current.cluster.ids <- c("B","DCs","Eosinophils","ILC_T","lowQ","Mono_Mac", "Neutrophils_0","Neutrophils_1",
                         "Neutrophils_2","Neutrophils_3", "Neutrophils_4","Neutrophils_5",
                         "Neutrophils_6")
new.cluster.ids <- c("B","DCs","Eosinophils","ILC_T","lowQ","Mono_Mac", "Neutrophils","Neutrophils",
                     "Neutrophils","Neutrophils", "Neutrophils","lowQ",
                     "lowQ")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "annotation")

### Mono/Mac 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Mono_Mac",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c("Mono_Mac_0","Mono_Mac_1","Mono_Mac_2","Mono_Mac_3","Mono_Mac_4",
                                        "Mono_Mac_5"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster", split.by = "condition")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC))

## marker genes from Kang B. et al 2020 Mucosal Immunology
#continuum: Ly6c high MHCII negative --> Ly6c high MHCII positive, Cx3cr1 intermediate, Il10, Tnfa, Il23 low --> 
#Ly6c int/low, MHCII high, Cx3cr1 high; cluster 1 expresses less Adgre1 compared to the others 
macs_gut_marker <- c("Adgre1" ,"Ly6c1","Ly6c2","Il10","Il23a", #summary of genes that change their expression 
                     "Thbs1", "Ccr2","Tgm2","Emilin2","Ifi27l2a", #Cluster 1 in Kang et al 
                     "Cx3cr1","Dst", "H2-M2","Dnmt3a","Vcam1", #Cluster 2 in Kang et al
                     "Apoe","Ms4a7","Cd63","Hpgds", #Cluster 3 in Kang et al (also Vcam1)
                     "Pgf","Mmp13","Dnase1l3","Acp5", "C1qb","C1qc", #Cluster 4 in Kang et al (also H2-M2)
                     "Ccl7","Ccl2","Pf4","Cd163", #Cluster 6 in Kang et al
                     "Hspa1a","Hes1","Hspa1b","Ccl4","Jun", #Cluster 7 in Kang et al 
                     "Rsad1","Ifit2","Isg15","Ifit1","Cxcl9" ,#Clsuter 11 in Kang et al
                     "Spp1","Fn1"
)

DotPlot_scCustom(sub_celltype, features = macs_gut_marker,dot.scale = 10, scale = TRUE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) +
  ggtitle("Kang et al marker") + theme(axis.text.x = element_text(angle = 90)) 

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

## rename
current.cluster.ids <- c("B","DCs","Eosinophils","ILC_T","lowQ","Mono_Mac_0","Mono_Mac_1","Mono_Mac_2","Mono_Mac_3","Mono_Mac_4",
                         "Mono_Mac_5","Neutrophils")
new.cluster.ids <- c("B","DCs","Eosinophils","ILC_T","lowQ","TAMs","lowQ","Monocytes","Macrophages","Macrophages",
                     "TRM","Neutrophils")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE, reduction = "umap.mnn")

### T/ILC 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "ILC_T",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.2)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "ILC_T_0","ILC_T_1","ILC_T_2","ILC_T_3","ILC_T_4"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

## marker genes from Smith P.M. et al 2011 Front. Microbiol and Seo G. et al 2020 Mucosal Immunology 
#ILC1 similar to Th1, ILC2 similar to Th2, ILC3 similar to Th17
#T-bet = Tbx21; Cd39 = Entpd1
T_marker <- c("Cd3g", #general marker of T cells 
                     "Il17a","Il22", "Entpd1",  #gamma delta T cells don't express CD4 and CD8, share markers with NK cells 
                     "Cd4","Foxp3", #Tregs (there could be two subsets, one expressing Il10 and one not)
                     "Cd8a", "Pdcd1", #CD8 T cells (also Tbx21 and Gzmb)
                     "Ccr7", #should only be high in naive CD8 T cells 
                     "Il21", #Th17 (also Il17a)
                     "Il10","Il13","Il5","Il4", #Th2
                     "Il12a", #Th1
                     "Mcpt1","Mcpt2", #Mast cells
                     "Tbx21", #ILC1 
                     "Gata3", #ILC2 (also Il4 and Il13 and Il5)
                     "Rorc",#ILC3 (also Il17, also Il22)
                     "Gzmb","Klrb1","Cd34", "Kit","Il15ra"#NK (also classified as killer ILCs) (also Cd8a)
)

DotPlot(sub_celltype, features = T_marker,dot.scale = 10, scale = TRUE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 15)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 15)) +
  ggtitle("") + theme(axis.text.x = element_text(angle = 90)) 

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

## rename
current.cluster.ids <- c("B", "DCs","Eosinophils","ILC_T_0","ILC_T_1","ILC_T_2","ILC_T_3","ILC_T_4","lowQ",
                         "Macrophages", "Monocytes", "Neutrophils","TAMs","TRM")
new.cluster.ids<- c("B", "DCs","Eosinophils","CD4_T","CD8_T","ILCs","mixed","lowQ","lowQ",
                    "Macrophages", "Monocytes", "Neutrophils","TAMs","TRM")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", reduction = "umap.mnn", label = TRUE)

### B cells 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "B",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, label = TRUE, group.by = "sub.cluster", reduction = "umap.mnn")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "B_0","B_1","B_2","B_3","B_4"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

#PD-L2 = Pdcd1lg2, Cd73 = Nt5e
b_marker <- c("Ighd","Cd19", #naive B cells
                     "Ighm", "Entpd1",  #mature B cells
                     "Cd80","Nt5e","Pdcd1lg2","Ighg1","Igha", #memory B cells 
              "Jchain","Slc3a2" #plasma cells, also Igha
)
DotPlot_scCustom(sub_celltype, features = b_marker,dot.scale = 10, scale = FALSE) + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 20)) +
  ggtitle("") + theme(axis.text.x = element_text(angle = 90)) 

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

## rename
current.cluster.ids <- c("B_0","B_1","B_2","B_3","B_4","CD4_T","CD8_T", "DCs","Eosinophils","ILCs", "lowQ",
                         "Macropahges", "Monocytes", "Neutrophils","TAMs","TRM")
new.cluster.ids <- c("mature_B","PCs","mixed","mixed","lowQ","CD4_T","CD8_T", "DCs","Eosinophils","ILCs", "lowQ",
                     "Macrophages", "Monocytes", "Neutrophils","TAMs","TRM")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", reduction = "umap.mnn", label = TRUE)

### remove lowQ and cycling cells 
Idents(subCl) <- "annotation"
obj <- subset(subCl, idents = c("CD4_T","CD8_T", "DCs","Eosinophils","ILCs",
                               "Macrophages","mature_B","Monocytes","Neutrophils","PCs","TAMs","TRM"))

DimPlot(obj, reduction = "umap.mnn", label = TRUE, group.by = "annotation",raster=FALSE)

##### Heatmap of annotation markers 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

markers <- c("Icos","Ctla4","Rora","Il2ra","Foxp3", #CD4 T
             "Trbc2","Ccl5","Ms4a4b","Nkg7","Cd8a",#CD8 T 
             "Hlf","Il17rb","Stab2","Ccl1","Calca",#ILCs
             "Scd1","Ighd","Pax5","Cd19","Cr2", #mature B 
             "Jchain","Igha","Iglv1","Ighg1","Tnfrsf17",#PCs
             "S100a8","S100a9","Lcn2","Slfn4","Cxcr2",#Neutrophils
             "Serpinb2","F5","Syne1","Ccr3","Alox15", #Eosinophils
             "Clec9a","Xcr1","Cd209a","Ccl22","Tbc1d4", #DCs
             "Ms4a6c","Ms4a4c","F10","Adora2b","Ifi205",#Monocytes
             "Vcam1","Ccl8","Gas6","Folr2","Igf1", #Macrophages
             "Arg1","Fn1","Spp1","Mmp12","Cxcl14",#TAMs
             "Pclaf","Birc5","Ccna2","Nuf2","Ube2c"#TRM
)
heatmap_goi_coi(obj, "annotation",markers,c("CD4_T","CD8_T","ILCs", "mature_B","PCs","Neutrophils","Eosinophils","DCs","Monocytes","Macrophages","TAMs","TRM"), 
                c(5,5,5,5,5,5,5,5,5,5,5,5),c("#270A7F",  "#26DFED","#6899C1", "#EDB20C","#DDED0C","#BD7FEA","#E22B17","#F20AB1",
                                             "#C0DBB4",  "#1E8209","#54EF0C","#D7EF92"),
                c(CD4_T="#270A7F", CD8_T= "#26DFED",ILCs="#6899C1", mature_B="#EDB20C",PCs="#DDED0C",Neutrophils="#BD7FEA",
                  Eosinophils="#E22B17",DCs="#F20AB1",Monocytes="#C0DBB4", Macrophages= "#1E8209",TAMs="#54EF0C", TRM= "#D7EF92"),F,T)

##### save object 
saveRDS(obj, "/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")
