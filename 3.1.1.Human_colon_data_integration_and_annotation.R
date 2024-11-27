########## This code integrates colon data from human biopsies and annotates clusters ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS(file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_combined.rds")

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
obj <- FindClusters(obj, resolution = 1, cluster.name = "mnn.clusters", algorithm = 2)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:15, reduction.name = "umap.mnn")
DimPlot(obj,reduction = "umap.mnn",group.by = "mnn.clusters",raster=FALSE)
obj <- JoinLayers(obj)

##### Annotation 

## DEG analysis
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
DimPlot(obj,reduction = "umap.mnn",raster=FALSE, group.by = "mnn.clusters",label = TRUE)
Idents(obj) <- "mnn.clusters"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

## heatmap of the top 5 DEGs 
markers <- c("ENSG00000286848","JCHAIN","IGHA2","IGKC","TNFRSF17", #cluster 0: PCs 
             "LILRB5","CD163L1","CD209","F13A1","IGSF21", #cluster 1 Macrophages 
             "IL7R","SPOCK2","ICOS","CD28","CCR4", #cluster 2: CD4+ T 
             "FFAR2","MARCHF3","DACH1","SYNE1","CLC", #cluster 3: Eosinophils
             "MME","PROK2","CXCR1","CMTM2","S100A12", #cluster 4: Neutrophils 
             "IGLC2","IGLC3","IGHA1","DERL3","IGLV6-57", #cluster 5: PCs
             "GZMA","GNLY","PRF1","GZMH","LINC02446", #cluster 6: CD8+ T 
             "FCN1","VCAN","CD300E","CLEC10A","CD1C", #cluster 7: DCs 
             "SPP1","APOC1","LPL","CHIT1","ENSG00000286618", #cluser 8: SPP1+ TAMs
             "CD24","CEACAM6","GPX2","TSPAN8","PIGR", #cluster 9: Epithelial cells 
             "CXCL8","IL1RN","HCAR3","HCAR2","KATNBL1", #cluter 10: Neutrophils 
             "ADAMDEC1","HAPLN1","STMN2","PCDH19","COL6A5", #cluster 11: Fibroblasts
             "IGHG1","IGHG2","IGHG3","IGHGP","IGHG4", #cluster 12: PCs
             "LEFTY1","C1orf210","PKD1-AS1","MUC5B","FGFBP1", #cluster 13: Epithelial cells 
             "RAPGEF2","RAB11FIP1","ANKRD28","SIK3","LINC-PINT", #cluster 14: Epithelial cells 
             "CPA3","TPSAB1","SIGLEC6","ADCYAP1","ENSG00000233968", #cluster 15: Mast cells 
             "OGN","THBS2","SFRP2","OMD","SFRP4", #cluster 16: Fibroblasts
             "ENSG00000290104" ,# cluster 17: PCs (the other from the top 5 DEGs are mathing the ones from PCs)
             "MS4A1","PAX5","CXCR5","FCRL1","NIBAN3", #cluster 18: mature B 
             "CD8A","ENSG00000290574","NKG7","CRTAM" ,#cluster 19: CD8+ T 
             "HSPA6","FER1L4","HSPA1L","ENSG00000285565","IFNG-AS1", #cluster 20: low quality 
             "CLDN5","SOX17","ACKR1","TM4SF18","SELE", #cluster 21: Endothelial 
             "CXCL3","CXCL5","IL6","TNFSF15","CSF3", #cluster 22: Macrophages 
             "LILRA4","VASH2","KRT5","LAMP3","TCL1A", #cluster 23: DCs
             "MS4A2","TPSB2","CALB2","TPSD1", #cluster 24: low quality 
             "MIR155HG","GNG7","TEX14","OSBPL3",  #cluster 25: low quality 
             "RHEX","CTSG",#clsuter 26: low quality 
             "MTUS2","ENSG00000225144", #cluster 27: low quality 
             "HSP90AA1","MALAT1","HSPA1A", #cluster 28: low qualtiy 
             "CEACAM8","LTF","CAMP","RETN","DEFA4", #cluster 29: Basophils 
             "CYGB","CCL8","CCL11","NTF3","WNT4",#cluster 30: low quality 
             "IGHM","BHLHA15","DTNB-AS1", #cluster 31: low quality 
             "MOXD1",#cluster 32: low quality 
             "ENSG00000240040","ANKRD36BP2"  , #cluster 33: low quality 
             "ALKAL2","PTGFR",  #clsuter 34: low quality 
             "CAVIN2","FCER1A","LYVE1","GGTA1" #cluster 35: low quality 
)

# heatmap per cluster and condition 
heatmap_goi_coi(obj, "mnn.clusters",markers,c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15",
                                             "c16","c17","c18","c19","cl20","c21","c22","c23","c24","c25","c26","c27","c28","c29",
                                             "c30","c31","c32","c33","c34","c35"), 
                c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,1,5,4,5,5,5,5,4,4,2,2,3,5,5,3,1,2,2,4),
                c("#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF",
                  "#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF",
                  "#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF",
                  "#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF","#EF0A3C",  "#0A41EF"),c
                (c0="#EF0A3C", c1= "#0A41EF",c2="#EF0A3C", c3= "#0A41EF",c4="#EF0A3C", c5= "#0A41EF",
                  c6="#EF0A3C", c7= "#0A41EF",c8="#EF0A3C", c9= "#0A41EF",c10="#EF0A3C", c11= "#0A41EF",
                  c12="#EF0A3C", c13= "#0A41EF", c14= "#EF0A3C",c15="#0A41EF", c16= "#EF0A3C",c17="#0A41EF", 
                  c18= "#EF0A3C",c19="#0A41EF",cl20="#EF0A3C",c21= "#0A41EF", c22="#EF0A3C", c23= "#0A41EF", 
                  c24= "#EF0A3C",c25="#0A41EF", c26= "#EF0A3C",c27="#0A41EF", c28= "#EF0A3C",c29="#0A41EF",
                  c30="#EF0A3C",c31= "#0A41EF",c32="#EF0A3C", c33= "#0A41EF", c34= "#EF0A3C",c35="#0A41EF"),F,F)

## heatmap of markers from Zheng X. et al 2022 Signal Trnasduc. Target. Ther. and Zilionis R. et al. 2019 Immunity 
markers <- c("MS4A2","KIT", #Mast cells
             "PECAM1","VWF", #Endothelial 
             "COL1A1","DCN", #Fibroblasts
             "FCGR3A","S100A12","LYZ", #Monocytes 
             "MSR1","MRC1","CD163","CD86","CD68", #M2 mac
             "IL1B","ITGAX", #M1 mac
             "LILRA4","IRF7","TCF4","IL3RA", #DCs
             "MZB1","IGLC2","IGHA1","IGHG3", #Plasma B 
             "MS4A1","CD79B","CD79A", #follicular B cells 
             "CD3G","CD3E","CD3D","PTPRC", #T cells 
             "KLRF1","KLRD1", #NKT
             "IL2RA","FOXP3", #Tregs
             "CD8B","CD8A", #CD8 T cells 
             "CXCR2","S100A8","CD14", #Neutrophils
             "CLC","CCR3","ADGRE1", #Eosinophils
             "HBB","HBA2","HBA1", #RBCs (red blood cells)
             "EPCAM", #Epithelial
             "CHGA",#Enteroendocrine
             "TFF3","MUC2", #Goblet
             "GUCA2B","GUCA2A","CA1","SLC26A3" #Enterocytes 
)

heatmap_goi_coi(obj, "mnn.clusters",markers,c("mast","endothelial","fibroblasts","monocytes","macM2","macM1","DCs", "plasmaB",
                                                "follicularB","Tcells","NKT","Treg","CD8_T","neutrophils","eosinophils",
                                                "RBCs", "epithelial","enterocrine","goblet","enterocytes"), 
                c(2,2,2,3,5,2,4,4,3,4,2,2,2,3,3,3,1,1,2,4),
                c("#EFA706",  "#EF9AEF","#E214E2",  "#309633","#1FE523",  "#6EF497","#025125","#D7DB14",  "#E2B44D",
                  "#0EA7DD",  "#0E76DD","#2EEDDF","#5353F2","#EF06A1","#EF062D", "#F4B3BD",  "#967D59","#EDD0A6","#72706E","#5B3B0C"),
                c(mast="#EFA706", endothelial= "#EF9AEF",fibroblasts="#E214E2",  monocytes= "#309633",macM2="#1FE523", macM1= "#6EF497",
                  DCs="#025125", plasmaB="#D7DB14", follicularB= "#E2B44D",Tcells="#0EA7DD", NKT= "#0E76DD",Treg="#2EEDDF", CD8_T= "#5353F2",
                  neutrophils="#EF06A1", eosinophils="#EF062D",RBCs="#F4B3BD", epithelial="#967D59", enterocrine= "#EDD0A6",goblet="#72706E",
                  enterocytes="#5B3B0C"),F,F)

## heatmap of markers from Salcher S. et al 2022 Cancer Cell 
markers <- c("BANK1","BLK", "CD19","CR2","CXCR5", #B cells
             "AOC1","HMSD","NCCRP1","TREML1", #DC mature
             "ACKR1","ADAMTSL1","ADCY4", #Endothleial cells
             "ADAMDEC1","ADORA3","AMDHD2","CCL13", #Macrophages 
             "ADCYAP1","CADPS", "CALB2","CDK15","CPA3", #Mast¨
             "KIR2DL1","KIR3DL1","LIM2", "MLC1", #NK
             "ADGRE3","ADGRG3","BTNL8","CMTM2","CXCR1","CXCR2","CYP4F3", #Neutrophils 
             "ABCB9","ALDH1L2","ALG14","ALG5", #Plasma cells 
             "ADAMTS12","ADAMTS2","ASPN", #Ctromal 
             "CD40LG","FXYD7","GTSCR1","IL2","IL7R", #CD4 T 
             "CCR5","CD8A","CD8B","CRTAM","GZMK", #CD8 T 
             "BATF","CCDC141","CCR4" ,"CCR8",#Treg
             "CLEC9A","XCR1", #cDC1
             "CD1C","CLEC10A", #cDC2
             "RAD54L"#pDCs
)

heatmap_goi_coi(obj, "mnn.clusters",markers,c("B","DC_mature","endothelial","macrophages","mast","NK","neutrophils", "PC",
                                                "stromal","CD4","CD8","Treg","cDC1","cDC2","pDCs"),c(5,4,3,4,5,4,7,4,3,5,5,4,2,2,1),
                c("#EFA706",  "#EF9AEF","#E214E2",  "#309633","#1FE523",  "#6EF497","#025125",
                  "#D7DB14",  "#E2B44D","#0EA7DD",  "#0E76DD","#2EEDDF", "#5353F2","#EF06A1","#EF062D"),
                c(B="#EFA706", DC_mature= "#EF9AEF",endothelial="#E214E2",  macrophages= "#309633",mast="#1FE523", NK= "#6EF497",
                  neutrophils="#025125", PC="#D7DB14", stromal= "#E2B44D",CD4="#0EA7DD", CD8= "#0E76DD",Treg="#2EEDDF", cDC1= "#5353F2",
                  cDC2="#EF06A1", pDCs="#EF062D"),F,F)

## nFeature and percent mito per cluster 
VlnPlot(obj, features = "nFeature_RNA",pt.size = 0)
VlnPlot(obj, features = "percent.mt",pt.size = 0)

## rename clusters 
current.cluster.ids <- c(0:35)
new.cluster.ids <- c("PC","Macs","T","Eosinophils","Neutrophils","PC","T",
                     "Macs","Macs","Epithelial","Neutrophils","Fibroblasts","PC",
                     "Epithelial","Epithelial","Mast","Fibroblasts","PC", "B_mature","T","lowQ", "Endothelial",
                     "Macs","Macs","lowQ","lowQ","lowQ","lowQ","lowQ","Basophils","lowQ","lowQ","lowQ","lowQ","lowQ","lowQ")
obj$annotation <- plyr::mapvalues(x = obj$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(obj, group.by = "annotation1", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of macrophages/DCs
Idents(obj) <- "annotation1"
subCl <- FindSubCluster(obj,cluster = "Macs",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.5)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Macs_0","Macs_1","Macs_2","Macs_3","Macs_4","Macs_5","Macs_6",
                                         "Macs_7","Macs_8","Macs_9","Macs_10"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

markers <- c("BATF3","XCR1", #cDC1
             "FCER1A","CD1C", #cDC2
             "LILRA4",#pDCs
             "CCR7","FSCN1", #activated DCs 
             "CD14","VCAN","S100A12", #Monocytes
             "LYZ","EREG","S100A6", #Mono-Mac transitory cells 
             "MRC1","C1QA","CD68","APOE","MARCO","NLRP3","CD163", #Macrophages
             "C1QC","SPP1" #SPP1+ TAMs
)
DotPlot(sub_celltype, features = markers,dot.scale = 10, scale = FALSE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

# rename
current.cluster.ids <- c("B_mature","Basophils","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Macs_0","Macs_1","Macs_2","Macs_3","Macs_4","Macs_5","Macs_6","Macs_7","Macs_8","Macs_9","Macs_10",
                         "Mast","Neutrophils","PC","T")
new.cluster.ids <- c("B_mature","Basophils","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Macrophages","TAMs","Monocytes","Macrophages","TAMs","Macrophages",
                     "Monocytes","DCs","lowQ","TAMs","DCs",
                     "Mast","Neutrophils","PC","T")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of T/ILCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "T",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.5)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "T_0","T_1","T_2","T_3","T_4","T_5","T_6","T_7"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- JoinLayers(sub_celltype)
sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

markers <- c("CD3G","CD4","RORC","IL7R", #Th17 cells, RORC also in ILC3, but no CD4 and CD3 there  
             "FOXP3", #T regs 
             "CD8A","CD8B",#CD8 T cell s
             "KLRG1","NMUR2" #ILC2
)
DotPlot(sub_celltype, features = markers,dot.scale = 10, scale = FALSE, assay = "RNA") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) + 
  theme(title = element_text(size = 20))+ theme(axis.text = element_text(size = 10)) + theme(axis.text.x = element_text(angle = 90)) 

# rename
current.cluster.ids <- c("B_mature","Basophils","DCs","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Mast","Monocytes","Neutrophils","PC",
                        "T_0","T_1","T_2","T_3","T_4","T_5","T_6","T_7","Macrophages","TAMs")
new.cluster.ids <- c("B_mature","Basophils","DCs","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Mast","Monocytes","Neutrophils","PC",
                    "CD4_T","CD4_T","CD8_T","CD8_T","lowQ","CD8_T","CD4_T","lowQ","Macrophages","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of B_mature
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "B_mature",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "B_mature_0","B_mature_1","B_mature_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- JoinLayers(sub_celltype)
sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

# rename
current.cluster.ids <- c("B_mature_0","B_mature_1","B_mature_2","Basophils","CD4_T","CD8_T",
                         "DCs","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts","lowQ","Mast","Monocytes","Neutrophils","PC",
                         "Macrophages","TAMs")
new.cluster.ids <- c("B_mature","low_Q","B_mature","Basophils","CD4_T","CD8_T",
                     "DCs","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts","lowQ","Mast","Monocytes","Neutrophils","PC",
                     "Macrophages","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of Fibroblasts 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Fibroblasts",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Fibroblasts_0","Fibroblasts_1","Fibroblasts_2","Fibroblasts_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- JoinLayers(sub_celltype)
sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

# rename
current.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                         "DCs","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts_0", "Fibroblasts_1", "Fibroblasts_2", "Fibroblasts_3","lowQ","Mast",
                         "Monocytes","Neutrophils","PC","Macrophages","TAMs")
new.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                     "DCs","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts", "Fibroblasts", "Fibroblasts", "lowQ","lowQ","Mast",
                     "Monocytes","Neutrophils","PC","Macrophages","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of Eosinophils
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Eosinophils",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Eosinophils_0","Eosinophils_1","Eosinophils_2"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- JoinLayers(sub_celltype)
sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

# rename
current.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                         "DCs","Endothelial","Eosinophils_0","Eosinophils_1","Eosinophils_2","Epithelial",
                         "Fibroblasts", "lowQ","low_Q","Mast",
                         "Monocytes","Neutrophils","PC","Macrophages","TAMs")
new.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                     "DCs","Endothelial","Eosinophils","Eosinophils","Eosinophils","Epithelial",
                     "Fibroblasts", "lowQ","lowQ","Mast",
                     "Monocytes","Neutrophils","PC","Macrophages","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of PCs
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "PC",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.1)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "PC_0","PC_1","PC_2","PC_3"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- JoinLayers(sub_celltype)
sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")

# rename
current.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                         "DCs","Endothelial","Eosinophils","Epithelial",
                         "Fibroblasts", "lowQ","Mast",
                         "Monocytes","Neutrophils","PC_0","PC_1","PC_2","PC_3","Macrophages","TAMs")
new.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                     "DCs","Endothelial","Eosinophils","Epithelial",
                     "Fibroblasts", "lowQ","Mast",
                     "Monocytes","Neutrophils","PC","PC","PC","lowQ","Macrophages","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### sub-clustering of Epithelial 
Idents(subCl) <- "annotation"
subCl <- FindSubCluster(subCl,cluster = "Epithelial",graph.name = "RNA_snn", 
                        subcluster.name = "sub.cluster",resolution = 0.3)
DimPlot(subCl, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

Idents(subCl) <- "sub.cluster"
sub_celltype <- subset(subCl,idents = c( "Epithelial_0","Epithelial_1","Epithelial_2",
                                         "Epithelial_3","Epithelial_4","Epithelial_5","Epithelial_6"))
DimPlot(sub_celltype, reduction = "umap.mnn", label = TRUE, group.by = "sub.cluster")

sub_celltype <- JoinLayers(sub_celltype)
sub_celltype <- NormalizeData(sub_celltype, normalization.method = "LogNormalize",
                              scale.factor = 10000,
                              margin = 1, assay = "RNA")
DefaultAssay(sub_celltype) <- "RNA"
markers <- FindAllMarkers(object = sub_celltype, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))

VlnPlot(sub_celltype, features = "nFeature_RNA")
VlnPlot(sub_celltype, features = "percent.mt")
FeaturePlot(sub_celltype, features = "MUC13", reduction = "umap.mnn")

# rename
current.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                     "DCs","Endothelial","Eosinophils","Epithelial_0","Epithelial_1","Epithelial_2","Epithelial_3",
                     "Epithelial_4","Epithelial_5","Epithelial_6",
                     "Fibroblasts", "lowQ","Mast",
                     "Monocytes","Neutrophils","PC","Macrophages","TAMs")
new.cluster.ids <- c("B_mature","Basophils","CD4_T","CD8_T",
                     "DCs","Endothelial","Eosinophils","Epithelial","lowQ","lowQ","Epithelial",
                     "lowQ","Epithelial","lowQ",
                     "Fibroblasts", "lowQ","Mast",
                     "Monocytes","Neutrophils","PC","Macrophages","TAMs")
subCl$annotation <- plyr::mapvalues(x = subCl$sub.cluster, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(subCl, group.by = "annotation", label = TRUE,raster=FALSE,reduction = "umap.mnn")

### remove lowQ cells 
Idents(subCl) <- "annotation"
obj <- subset(subCl, idents = c("B_mature","Basophils","CD4_T","CD8_T",
                               "DCs","Endothelial","Eosinophils","Epithelial",
                               "Fibroblasts","Mast","Monocytes","Neutrophils","PC","Macrophages","TAMs"))

##### Heatmap of annotation markers 
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(obj) <- "annotation"
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
write.csv(markers, "/scratch/khandl/eos_human/data_files/DEGs/DEGs_annotated_clusters_patient_biopsies.csv")
View(markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC))

markers <- c("IL7R","ICOS","SPOCK2","CCR4","KLRB1", #CD4 T
             "GZMA","CD8A","NKG7","GZMB","GZMK",#CD8 T 
             "MS4A1","BANK1","CXCR5","CD19","CR2", #mature B 
             "JCHAIN","IGHA2","IGHA1","IGLC3","IGHG1",#PCs
             "BPI","CEACAM8","LTF","DEFA3","DEFA4", #Basophils
             "FCGR3B","CXCR2","S100A8","S100A9","HCAR3",#Neutrophils
             "FFAR2","SYNE1","CLC","DACH1","MARCHF3", #Eosinophils
             "CPA3","TPSAB1","MS4A2","TPSB2","SIGLEC6",#Mast
             "CD1C","CLEC10A","FCER1A","PKIB","FLT3", #DCs
             "FCN1","VCAN","EREG","THBS1","CD300E","FCN1", #Monocytes
             "CD209","F13A1","LILRB5","DNASE1L3","SELENOP","C1QC","C1QB", #Macrophages
             "SPP1","APOC1","TREM2","LPL","MARCO"#TAMs
)
heatmap_goi_coi(obj, "annotation",markers,c("CD4_T","CD8_T","mature_B","PCs","Basophils","Neutrophils","Eosinophils","Mast","DCs","Monocytes","Macrophages","TAMs"), 
                c(5,5,5,5,5,5,5,5,5,6,7,5),c("#270A7F",  "#26DFED","#EDB20C","#DDED0C",  "#7C7C79","#BD7FEA","#E22B17",  "#56544F","#F20AB1",
                                             "#C0DBB4",  "#1E8209","#54EF0C"),
                c(CD4_T="#270A7F", CD8_T= "#26DFED",mature_B="#EDB20C",PCs="#DDED0C", Basophils= "#7C7C79",Neutrophils="#BD7FEA",
                  Eosinophils="#E22B17", Mast= "#56544F",DCs="#F20AB1",Monocytes="#C0DBB4", Macrophages= "#1E8209",TAMs="#54EF0C"),F,T)

##### save object 
saveRDS(obj, "/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")


