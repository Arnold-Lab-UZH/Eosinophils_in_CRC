########## This code analyses differences in cell type proportions between tumor and control (NAT) ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")

#only consider CD45 positive cells 
Idents(obj) <- "annotation"
obj <- subset(obj, idents = c("B_mature","Basophils","CD4_T","CD8_T","DCs","Eosinophils","Mast","Monocytes",
                                   "Macrophages","Neutrophils","PC","TAMs"))

##### plot umap split by condition 
p <- DimPlot(obj, group.by = "annotation", label = FALSE,
             cols = c("#EDB20C","#7C7C79","#270A7F", "#26DFED","#F20AB1","#E22B17"
                      ,"#1E8209","#56544F", "#C0DBB4","#BD7FEA","#DDED0C","#54EF0C"),
             raster=FALSE, reduction = "umap.mnn", split.by = "tissue")
ggsave(file = "/scratch/khandl/eos_human/figures/anno/umap_anno.svg", width = 8, height = 6, plot = p)

##### plot proportions in barplot 
create_table_cell_type_prop(obj, "tissue","annotation","/scratch/khandl/eos_human/figures/","TME")
df <- read.csv("/scratch/khandl/eos_human/figures/TME_proportions_tissue_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:13)) 

#define order 
df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c( "Eosinophils","B_mature","PC","Basophils","Neutrophils",
                                                                    "CD4_T","CD8_T","DCs","Mast","Monocytes","Macrophages",
                                                                    "TAMs" )))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E22B17","#EDB20C","#DDED0C","#7C7C79", "#BD7FEA","#270A7F","#26DFED","#F20AB1","#56544F",
                               "#C0DBB4", "#1E8209", "#54EF0C")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/figures/anno/TME_hs_proportions.svg", width = 12, height = 6, plot = p)

##### barplot per patient
create_table_cell_type_prop(obj, "condition","annotation","/scratch/khandl/eos_human/figures/","TME")
df <- read.csv("/scratch/khandl/eos_human/figures/TME_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:13)) 

#define order 
df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c( "Eosinophils","B_mature","PC","Basophils","Neutrophils",
                                                                    "CD4_T","CD8_T","DCs","Mast","Monocytes","Macrophages",
                                                                    "TAMs" )))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E22B17","#EDB20C","#DDED0C","#7C7C79", "#BD7FEA","#270A7F","#26DFED","#F20AB1","#56544F",
                               "#C0DBB4", "#1E8209", "#54EF0C")) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_human/figures/anno/TME_hs_proportions_per_patient.svg", width = 12, height = 6, plot = p)

##### statistical testing  
#first control, then the condition
cell_type_prop_stats(obj,"annotation","tissue_ctrl","tumor","tissue",1.41,
                     "/scratch/khandl/eos_human/figures/cell_type_prop/tissue_ctrl_vs_tumor_stat.svg") 

##### Barplot to compare eos, neutrophils and TAMs 
df <- read.csv("/scratch/khandl/eos_human/figures/TME_proportions_tissue_annotation.csv", header = TRUE)

barplot_cell_type_oi_2cond(df,c("X","Neutrophils"),c("#BCA6C9","#690E9E"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/human_barplot_neut.svg")
barplot_cell_type_oi_2cond(df,c("X","Eosinophils"),c("#E08282","#930F0F"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/human_barplot_eos.svg")
barplot_cell_type_oi_2cond(df,c("X","TAMs"),c("#C5EFA3","#4A9B0A"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/human_barplot_tams.svg")
barplot_cell_type_oi_2cond(df,c("X","PC"),c("#E0E0A2","#D8D805"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/_human_barplot_PC.svg")
barplot_cell_type_oi_2cond(df,c("X","B_mature"),c("#EACD8C","#A5750A"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/human_barplot_matureB.svg")
