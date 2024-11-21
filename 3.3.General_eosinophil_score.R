########## This code applies a general eosinophil score based on DEGs from mouse data to human data ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load data 
### Reference data (mouse data)
blood_mouse <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_blood_annotated.rds")
colon_mouse <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")

### Test datasets (human) data 
blood_human <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_blood_annotated.rds")
colon_human <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")

##### colon score 
### extract the top 50 DEGs from the mouse eosinophil cluster 
mouse_colon_DEGs <- NormalizeData(colon_mouse, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(colon_mouse) <- "annotation"
markers <- FindAllMarkers(object = colon_mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
markers_eos <- markers[markers$cluster %in% "Eosinophils",]
markers_top50 <- markers_eos %>% top_n(n =50, wt = avg_log2FC)

### convert to human symbols
rownames(markers_top50) <- markers_top50$gene
eos_markers_colon_human_orthologues <- convert_orthologs(gene_df = markers_top50,
                                                         gene_input = "rownames", 
                                                         gene_output = "rownames", 
                                                         input_species = "mouse",
                                                         output_species = "human",
                                                         non121_strategy = "drop_both_species",
                                                         method = "gprofiler") 
### add the score 
colon_human <-AddModuleScore(colon_human, features= list(rownames(eos_markers_colon_human_orthologues)),name = "Mouse_eos_score")

### plot feature plot
p <- FeaturePlot(object = colon_human, features = "Mouse_eos_score1", pt.size = 1, reduction = "umap.mnn") + theme(legend.position = "right") + 
  scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,1.2))
ggsave(file = "/scratch/khandl/eos_human/figures/eos_score/human_tme_eos_score.pdf", width = 6, height = 5, plot = p)

### plot vlnplot
p <- VlnPlot(colon_human, features= "Mouse_eos_score1", group.by = "annotation", pt.size = 0, 
             c("#EDB20C","#7C7C79","#270A7F", "#26DFED","#F20AB1","#BCB2A2","#E22B17","#AA8C50",
               "#422D03","#1E8209","#56544F", "#C0DBB4","#BD7FEA","#DDED0C","#54EF0C"))+  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave("/scratch/khandl/eos_human/figures/eos_score/colon_score_vln.svg", width = 8, height = 8, plot = p)

### statistics 
Idents(colon_human) <- "annotation"
eos <- subset(colon_human, idents = "Eosinophils")
baso <- subset(colon_human, idents = "Basophils")
p_val <- wilcox.test(eos$Mouse_eos_score1, baso$Mouse_eos_score1, alternative = "two.sided")
print(p_val)

##### blood score 
### extract the top 50 DEGs from mouse eosinophil clusters 
mouse_blood_DEGs <- NormalizeData(blood_mouse, normalization.method = "LogNormalize",scale.factor = 10000,margin = 1, assay = "RNA")
Idents(blood_mouse) <- "annotation"
markers <- FindAllMarkers(object = blood_mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
markers_eos <- markers[markers$cluster %in% "Eosinophils",]
markers_top50 <- markers_eos %>% top_n(n =50, wt = avg_log2FC)

### convert to human symbols
rownames(markers_top50) <- markers_top50$gene
eos_markers_blood_human_orthologues <- convert_orthologs(gene_df = markers_top50,
                                                         gene_input = "rownames", 
                                                         gene_output = "rownames", 
                                                         input_species = "mouse",
                                                         output_species = "human",
                                                         non121_strategy = "drop_both_species",
                                                         method = "gprofiler") 

### add the score 
blood_human <-AddModuleScore(blood_human, features= list(rownames(eos_markers_blood_human_orthologues)),name = "Mouse_eos_score")

### plot feature plot
p <- FeaturePlot(object = blood_human, features = "Mouse_eos_score1", pt.size = 1, reduction = "umap.mnn") + theme(legend.position = "right") + 
  scale_color_gradientn( colours = c('grey', 'darkred'),  limits = c(0,1.2))
ggsave(file = "/scratch/khandl/eos_human/figures/eos_score/human_blood_eos_score.pdf", width = 6, height = 5, plot = p)

### plot vlnplot
p <- VlnPlot(blood_human, features= "Mouse_eos_score1", group.by = "annotation", pt.size = 0, 
             c("#E22B17", "#BD7FEA", "#54EF0C","#4166EF")) +  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
ggsave("/scratch/khandl/eos_human/figures/eos_score/blood_score_vln.svg", width = 8, height = 8, plot = p)

### statistics
Idents(blood_human) <- "annotation"
eos <- subset(blood_human, idents = "Eosinophils")
neut <- subset(blood_human, idents = "Neutrophils")
p_val <- wilcox.test(eos$Mouse_eos_score1, neut$Mouse_eos_score1, alternative = "two.sided")
print(p_val)

