########## This code annotates human eosinophil clusters after integration ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
human_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated.rds")
mouse_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_eosinophils_integrated_annotated.rds")

#only take wild type samples 
Idents(mouse_eos) <- "condition"
mouse_eos_wt <- subset(mouse_eos, idents = c("adult_colon_wt","tumor_wt","blood_wt"))

##### identify the top 50 DEGs from eosinoophil subtypes to apply a module score approach 
## extract top 50 DEGs per subtype 
mouse_eos_wt <- NormalizeData(mouse_eos_wt , normalization.method = "LogNormalize", scale.factor = 10000,margin = 1, assay = "RNA")
Idents(mouse_eos_wt) <- "annotation"
markers <- FindAllMarkers(object = mouse_eos_wt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA", slot = "data")
markers_df <- as.data.frame(markers %>% group_by(cluster) %>% top_n(n =50, wt = avg_log2FC))

## convert the top 50 genes to human symbols 
#A eos
a_eos_markers <- markers_df[markers_df$cluster %in% "A_eos",]
rownames(a_eos_markers) <- a_eos_markers$gene
a_eos_human_orthologues <- convert_orthologs(gene_df = a_eos_markers,
                                             gene_input = "rownames", 
                                             gene_output = "rownames", 
                                             input_species = "mouse",
                                             output_species = "human",
                                             non121_strategy = "drop_both_species",
                                             method = "gprofiler") 
a_eos_human_orthologues <- rownames(a_eos_human_orthologues)

# basal eos
b_eos_markers <- markers_df[markers_df$cluster %in% "B_eos",]
rownames(b_eos_markers) <- b_eos_markers$gene
b_eos_human_orthologues <- convert_orthologs(gene_df = b_eos_markers,
                                             gene_input = "rownames", 
                                             gene_output = "rownames", 
                                             input_species = "mouse",
                                             output_species = "human",
                                             non121_strategy = "drop_both_species",
                                             method = "gprofiler") 
b_eos_human_orthologues <- rownames(b_eos_human_orthologues)

## add the score to the seurat object 
human_eos <-AddModuleScore(human_eos, features= list(a_eos_human_orthologues),name = "A_eos_score")
human_eos <-AddModuleScore(human_eos, features= list(b_eos_human_orthologues),name = "B_eos_score")

### plot module scores by mnn clusters 
Idents(human_eos) <- "mnn.clusters"
p <- VlnPlot(human_eos, features= "A_eos_score1", group.by = "mnn.clusters", pt.size = 0, cols =  c("#CA13EA","#EAA10F")) +  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") 
ggsave(file = "/scratch/khandl/eos_human/eos_anno/A_vln.svg", width = 8, height = 6, plot = p)

p <- VlnPlot(human_eos, features= "B_eos_score1", group.by = "mnn.clusters", pt.size = 0, cols =  c("#CA13EA","#EAA10F")) +  
  theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") 
ggsave(file = "/scratch/khandl/eos_human/eos_anno/B_vln.svg", width = 8, height = 6, plot = p)

#### statistical test 
# A-eos
Idents(human_eos) <- "mnn.clusters"
first <- subset(human_eos, idents = 1)
second <- subset(human_eos, idents = 2)
p_val <- wilcox.test(first$A_eos_score1, second$A_eos_score1, alternative = "two.sided")
p_val$p.value

# B-eos
Idents(human_eos) <- "mnn.clusters"
first <- subset(human_eos, idents = 1)
second <- subset(human_eos, idents = 2)
p_val <- wilcox.test(first$A_eos_score1, second$A_eos_score1, alternative = "two.sided")
p_val$p.value

#### rename clusters 
Idents(human_eos) <- "mnn.clusters"
current.cluster.ids <- c(1:2)
new.cluster.ids <- c("B_eos","A_eos")
human_eos$annotation <- plyr::mapvalues(x = human_eos$mnn.clusters, from = current.cluster.ids, to = new.cluster.ids)
p <- DimPlot(human_eos,reduction = "umap.mnn",group.by = "annotation",raster=FALSE, label = TRUE, cols = c("#10A069","#E81818"))
ggsave(file = "/scratch/khandl/eos_human/eos_anno/umap_annotated.svg", width = 8, height = 6, plot = p)

##### proportion plot per condition 
## per tissue 
create_table_cell_type_prop(human_eos, "tissue","annotation","/scratch/khandl/eos_human/eos_anno/","obj")
df <- read.csv("/scratch/khandl/eos_human/eos_anno/obj_proportions_tissue_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E81818","#10A069")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/eos_anno/cell_type_prop.svg", width = 12, height = 6, plot = p)

## per patient and condition 
Idents(human_eos) <- "tissue"
sub <- subset(human_eos, idents = "tumor") 
create_table_cell_type_prop(sub, "condition","annotation","/scratch/khandl/eos_human/eos_anno/","tumor")
df <- read.csv("/scratch/khandl/eos_human/eos_anno/tumor_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E81818","#10A069")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/eos_anno/cell_type_prop_per_patinet_tumor.svg", width = 12, height = 6, plot = p)

Idents(human_eos) <- "tissue"
sub <- subset(human_eos, idents = "tissue_ctrl") 
create_table_cell_type_prop(sub, "condition","annotation","/scratch/khandl/eos_human/eos_anno/","tissue_ctrl")
df <- read.csv("/scratch/khandl/eos_human/eos_anno/tissue_ctrl_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E81818","#10A069")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/eos_anno/cell_type_prop_per_patinet_tissue_ctrl.svg", width = 12, height = 6, plot = p)

Idents(human_eos) <- "tissue"
sub <- subset(human_eos, idents = "blood_healhty") 
create_table_cell_type_prop(sub, "condition","annotation","/scratch/khandl/eos_human/eos_anno/","blood_healthy")
df <- read.csv("/scratch/khandl/eos_human/eos_anno/blood_healthy_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:3)) 

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E81818","#10A069")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/eos_anno/cell_type_prop_per_patinet_bloodH.svg", width = 12, height = 6, plot = p)

##### statistics between tissue control and tumor
#first control, then the condition
Idents(human_eos) <- "annotation"
human_eos <- subset(human_eos, idents = c("B_eos","A_eos"))
cell_type_prop_stats(human_eos,"annotation","tissue_ctrl","tumor","tissue",1.41,
                     "/scratch/khandl/eos_human/eos_anno/stats_tissue_ctrl_vs_tumor.svg") 

##### save object 
saveRDS(human_eos,"/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated_annotated.rds")
