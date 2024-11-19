########## This code integrated mouse eosinophils from tumor, control and blood ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_CRC/1.Packages_and_functions.R")

##### read in R objects 
colon <- readRDS(file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")
blood <- readRDS(file =  "/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_blood_annotated.rds")

# extract eosinophils 
Idents(colon) <- "annotation"
colon_eos <- subset(colon, idents = c("Eosinophils"))

Idents(blood) <- "annotation"
blood_eos <- subset(blood, idents = c("Eosinophils"))

# only extract wt samples and exclude disseminated --> only 7 cells 
Idents(colon_eos) <- "condition"
colon_eos <- subset(colon_eos, idents = c("adult_colon_wt","tumor_wt"))

# merge datasets 
obj <- merge(colon_eos, blood_eos)
obj <- JoinLayers(obj)

### load reference 
adult_eosSS <- readRDS("/data/khandl/common/Nature_paper_data/eosinophils_steadystate.rds")
adult_eosSS$annotation <- adult_eosSS$seurat_clusters

#remove spleen, stomach and small intestine, keep bonemarrow for circ eos  
Idents(adult_eosSS) <- "orig.ident"
adult_eosSS <- subset(adult_eosSS, idents = c("bonemarrow","blood","colon"))

##### Pre-process
adult_eosSS <- NormalizeData(adult_eosSS)
adult_eosSS <- FindVariableFeatures(adult_eosSS)
adult_eosSS <- ScaleData(adult_eosSS,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
adult_eosSS <- RunPCA(adult_eosSS)
adult_eosSS <- FindNeighbors(adult_eosSS, reduction = "pca", dims = 1:30)
adult_eosSS <- FindClusters(adult_eosSS, resolution = 2, algorithm = 2)
adult_eosSS <- RunUMAP(adult_eosSS, reduction = "pca", dims = 1:30,return.model=TRUE)
adult_eosSS$condition <- adult_eosSS$orig.ident
DimPlot(adult_eosSS,reduction = "umap",group.by = "annotation",combine = FALSE, label.size = 2)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5, algorithm = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30,return.model=TRUE)
DimPlot(obj,reduction = "umap",combine = FALSE, label.size = 2, split.by = "condition")

##### Label transfer 
query <- obj
query <- NormalizeData(query)
anchors <- FindTransferAnchors(reference = adult_eosSS, query = query, dims = 1:30,reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = adult_eosSS$annotation, dims = 1:30)
query <- AddMetaData(query, metadata = predictions)
#MapQuery is a combined function for TransferData(), IntegrateEmbeddings() and ProjectUMAP()
query <- MapQuery(anchorset = anchors, reference = adult_eosSS, query = query,
                  refdata = list(celltype = "annotation"), reference.reduction = "pca", reduction.model = "umap")
p1 <- DimPlot(adult_eosSS, reduction = "umap", group.by = "annotation", label = TRUE, label.size = 3,pt.size = 0.5,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,pt.size = 0.5,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

# plot prediction score per condition
p <- VlnPlot(query, features = "predicted.celltype.score", group.by = "condition")
ggsave(file = "/scratch/khandl/eos_tumor/eosinophils/prediction_score.svg", plot = p, width = 10, height = 6)
#query$celltype <- ifelse(query$prediction.score.max > 0.5, query$predicted.id, "non.defined")

#add a meta.data column combining annotation and predicted.celltype 
adult_eosSS$annotation <- adult_eosSS$annotation
query$annotation <- query$predicted.celltype
#merge objects with umap and ref.umap based on the adult_ss 
integrated <- merge(adult_eosSS, query)
integrated[["umap"]] <- merge(adult_eosSS[["umap"]], query[["ref.umap"]])
DimPlot(integrated, group.by = "condition")

## rename non-defined A-eos 
current.cluster.ids <- c("basal eosinophils", "circulating eosinophils","eosinophil progenitors","immature eosinophils",
                         "intestinal eosinophils")
new.cluster.ids <- c("B_eos", "circ_eos","precursor_eos","immature_eos",
                     "A_eos")
integrated$annotation <- plyr::mapvalues(x = integrated$annotation, from = current.cluster.ids, to = new.cluster.ids)

# Plot tumor cells cells as an overlay to the reference and the reference separate
Idents(integrated) <- "condition"
sub <- subset(integrated, idents = c("blood","bonemarrow","colon"))
p <- DimPlot(sub, group.by = "annotation1",cols = c("#E81818","#10A069","#E8E81A", "#E88A1A","#26DFED"), pt.size = 1)
ggsave(file = "/scratch/khandl/eos_tumor/eosinophils/umap_annotated.svg", plot = p, width = 10, height = 6)

#highlihght tumor cells with larger dots 
Idents(integrated) <- "condition"
sub <- subset(integrated, idents = c("blood","bonemarrow","colon", "tumor_wt"))
tumor_cells <- WhichCells(sub, idents = c("tumor_wt"))
p <- DimPlot(sub, label=T, group.by="condition", cells.highlight= list(tumor_cells), cols.highlight = c( "#EA0A30"), cols= "#A39F9F")
ggsave(file = "/scratch/khandl/eos_tumor/eosinophils/highlight_tumor_cell_highlihght.svg", plot = p, width = 10, height = 6)

##### compare proportions between conditions 
create_table_cell_type_prop(integrated, "condition","annotation","/scratch/khandl/eos_tumor/eosinophils/","integrated")
df <- read.csv("/scratch/khandl/eos_tumor/eosinophils/integrated_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:6)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E81818","#10A069","#E8E81A", "#E88A1A","#26DFED")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_tumor/eosinophils/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### QC per condition 
Idents(integrated) <- "condition"
obj <- subset(integrated, idents = c("adult_colon_wt","blood_wt","tumor_wt"))
p <- VlnPlot(obj, features = "nFeature_RNA", pt.size = 0)
ggsave(file = "/scratch/khandl/eos_tumor/QC/nFeature_mouse_anno.svg", width = 8, height = 6, plot = p)
p <- VlnPlot(obj, features = "nCount_RNA", pt.size = 0)
ggsave(file = "/scratch/khandl/eos_tumor/QC/nCoung_mouse_anno.svg", width = 8, height = 6, plot = p)
p <- VlnPlot(obj, features = "percent.mt", pt.size = 0)
ggsave(file = "/scratch/khandl/eos_tumor/QC/percent_mt_mouse_anno.svg", width = 8, height = 6, plot = p)

### save R object 
saveRDS(integrated,"/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_eosinophils_integrated_annotated.rds")
