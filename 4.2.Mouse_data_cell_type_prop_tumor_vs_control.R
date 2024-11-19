########## This code analyses differences in cell type proportions between colon samples (tumor vs. control) ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")

##### plot umap split by condition 
p <- DimPlot(obj, group.by = "annotation", pt.size = 0.1, label = FALSE, label.size = 3,reduction = "umap.mnn",
             cols = c("#270A7F", "#26DFED","#F20AB1",
                      "#E22B17","#6899C1","#1E8209","#EDB20C", "#C0DBB4","#BD7FEA", "#DDED0C",
                      "#54EF0C","#D7EF92"), split.by = "condition")
ggsave("/scratch/khandl/eos_tumor/figures/anno/umap_annotated.svg", width = 8, height = 5, plot = p)

##### plot proportions in barplot 
create_table_cell_type_prop(obj, "condition","annotation","/scratch/khandl/eos_tumor/data_files/","tumor_immune")
df <- read.csv("/scratch/khandl/eos_tumor/data_files/tumor_immune_proportions_condition_annotation.csv", header = TRUE)

df_plotting <- create_table_cell_type_prop_table_for_plot(df,c(2:13)) 

#define order
df_plotting <- within(df_plotting, cell_types <- factor(cell_types, 
                                                        levels = c("Eosinophils","mature_B","PCs","Neutrophils",
                                                                   "CD4_T","CD8_T","ILCs" ,"DCs","Monocytes","Macrophages","TRM",
                                                                   "TAMs")))

p <- ggplot(data=df_plotting, aes(x=sample, y=proportion, fill = cell_types)) +
  geom_bar(stat="identity", position = "fill" ) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.title= element_text(size = 25)) + theme(legend.title = element_text(size = 30), legend.text = element_text(size = 30)) + 
  xlab("Sample") + ylab("Cell type proportion") + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values=  c("#E22B17", "#EDB20C","#DDED0C","#BD7FEA", "#270A7F","#26DFED","#6899C0","#F281CC",
                               "#C0DBB4","#1E8209","#D7EF92", "#54EF0C"
  )) + coord_flip() + 
  theme_classic(base_size = 25) 
ggsave("/scratch/khandl/eos_tumor/figures/anno/cell_type_prop.svg", width = 12, height = 6, plot = p)

##### statistical testing  
#first control, then the condition
cell_type_prop_stats(obj,"annotation","tumor_wt","tumor_phil","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_tumor_wt_vs_tumor_phil.svg") 

cell_type_prop_stats(obj,"annotation","disseminated_wt","disseminated_phil","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_mets_wt_vs_mets_phil.svg") 

cell_type_prop_stats(obj,"annotation","adult_colon_wt","adult_colon_phil","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_adult_colon_wt_vs_adult_colon_phil.svg") 

cell_type_prop_stats(obj,"annotation","adult_colon_wt","tumor_wt","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_adult_colon_wt_vs_tumor_wt.svg") 

cell_type_prop_stats(obj,"annotation","disseminated_wt","tumor_wt","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_mets_wt_tumor_wt.svg") 

cell_type_prop_stats(obj,"annotation","adult_colon_phil","tumor_phil","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_adult_colon_phil_vs_tumor_phil.svg") 

cell_type_prop_stats(obj,"annotation","disseminated_phil","tumor_phil","condition",1.41,
                     "/scratch/khandl/eos_tumor/figures/cell_type_prop/stats_mets_phil_tumor_phil.svg") 

##### Barplot to compare eos, neutrophils and TAMs  control vs. tumor 
df <- read.csv("/scratch/khandl/eos_tumor/data_files/tumor_immune_proportions_condition_annotation.csv", header = TRUE)
df2 <- df[df$X %in% c("adult_colon_wt","tumor_wt"),]

barplot_cell_type_oi_2cond(df2,c("X","Neutrophils"),c("#BCA6C9","#690E9E"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/barplot_neut.svg")
barplot_cell_type_oi_2cond(df2,c("X","Eosinophils"),c("#E08282","#930F0F"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/barplot_eos.svg")
barplot_cell_type_oi_2cond(df2,c("X","TAMs"),c("#C5EFA3","#4A9B0A"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/barplot_tams.svg")
barplot_cell_type_oi_2cond(df2,c("X","PCs"),c("#E0E0A2","#D8D805"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/barplot_PC.svg")
barplot_cell_type_oi_2cond(df2,c("X","mature_B"),c("#EACD8C","#A5750A"),"/scratch/khandl/eos_tumor/figures/cell_type_prop/barplot_matureB.svg")

##### barplot to compare WT vs. PHIL in control, tumor and disseminated 
## eosinophils 
df <- read.csv("/scratch/khandl/eos_tumor/data_files/tumor_immune_proportions_condition_annotation.csv", header = TRUE)
df <- df[df$X %in% c("adult_colon_wt", "adult_colon_phil","tumor_wt","tumor_phil","disseminated_wt","disseminated_phil"),colnames(df) %in% c("X","Eosinophils")]
df1 <- df[df$X %in% c("adult_colon_wt","tumor_wt","disseminated_wt"),]
df1$condition <- "wt"
df2 <- df[df$X %in% c("adult_colon_phil","tumor_phil","disseminated_phil"),]
df2$condition <- "phil"

df <- rbind(df1, df2)

# Grouped
p <- ggplot(df, aes(fill=condition, y=Eosinophils, x=X)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 10)) + 
  theme(axis.title= element_text(size = 10))  + 
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))  + 
  xlab("Sample") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#E08082","#951C1F"))
print(p)
ggsave("/scratch/khandl/eos_tumor/grouped_barplot_eos.svg", width = 4, height = 5, plot = p)

## TAMs 
df <- read.csv("/scratch/khandl/eos_tumor/data_files/tumor_immune_proportions_condition_annotation.csv", header = TRUE)
df <- df[df$X %in% c("adult_colon_wt", "adult_colon_phil","tumor_wt","tumor_phil","disseminated_wt","disseminated_phil"),colnames(df) %in% c("X","TAMs")]
df1 <- df[df$X %in% c("adult_colon_wt","tumor_wt","disseminated_wt"),]
df1$condition <- "wt"
df2 <- df[df$X %in% c("adult_colon_phil","tumor_phil","disseminated_phil"),]
df2$condition <- "phil"

df <- rbind(df1, df2)

# Grouped
p <- ggplot(df, aes(fill=condition, y=TAMs, x=X)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 10)) + 
  theme(axis.title= element_text(size = 10))  + 
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))  + 
  xlab("Sample") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#C7E1A3","#499C45"))
print(p)
ggsave("/scratch/khandl/eos_tumor/grouped_barplot_tams.svg", width = 4, height = 5, plot = p)

## Macrophages 
df <- read.csv("/scratch/khandl/eos_tumor/data_files/tumor_immune_proportions_condition_annotation.csv", header = TRUE)
df <- df[df$X %in% c("adult_colon_wt", "adult_colon_phil","tumor_wt","tumor_phil","disseminated_wt","disseminated_phil"),colnames(df) %in% c("X","Macrophages")]
df1 <- df[df$X %in% c("adult_colon_wt","tumor_wt","disseminated_wt"),]
df1$condition <- "wt"
df2 <- df[df$X %in% c("adult_colon_phil","tumor_phil","disseminated_phil"),]
df2$condition <- "phil"

df <- rbind(df1, df2)

# Grouped
p <- ggplot(df, aes(fill=condition, y=Macrophages, x=X)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 10)) + 
  theme(axis.title= element_text(size = 10))  + 
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))  + 
  xlab("Sample") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#57715E","#116C37"))
print(p)
ggsave("/scratch/khandl/eos_tumor/grouped_barplot_macs.svg", width = 4, height = 5, plot = p)

## Monocytes 
df <- read.csv("/scratch/khandl/eos_tumor/data_files/tumor_immune_proportions_condition_annotation.csv", header = TRUE)
df <- df[df$X %in% c("adult_colon_wt", "adult_colon_phil","tumor_wt","tumor_phil","disseminated_wt","disseminated_phil"),colnames(df) %in% c("X","Monocytes")]
df1 <- df[df$X %in% c("adult_colon_wt","tumor_wt","disseminated_wt"),]
df1$condition <- "wt"
df2 <- df[df$X %in% c("adult_colon_phil","tumor_phil","disseminated_phil"),]
df2$condition <- "phil"

df <- rbind(df1, df2)

# Grouped
p <- ggplot(df, aes(fill=condition, y=Monocytes, x=X)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +
  theme(axis.text = element_text(size = 10)) + 
  theme(axis.title= element_text(size = 10))  + 
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))  + 
  xlab("Sample") + ylab("Cell type proportion") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=  c("#C6CCC2","#C2DBB4"))
print(p)
ggsave("/scratch/khandl/eos_tumor/grouped_barplot_mono.svg", width = 4, height = 5, plot = p)









