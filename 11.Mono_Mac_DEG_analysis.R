########## This code analyses DEGs betweeen tumor and control of the mono/mac populations (human and mouse)  ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
human <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")
mouse <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")

##### run DEG analysis 
## cell types to test 
celltypes <- c("TAMs","Macrophages","Monocytes")

### human data tumor vs. tissue ctrl 
for(i in celltypes) {
  Idents(human) <- "annotation"
  sub <- subset(human, idents = i)
  Idents(sub) <- "tissue"
  DEG_to_csv_two_cond(sub,"tumor","tissue_ctrl",FALSE,0.25,
                      paste0("/scratch/khandl/eos_human/data_files/DEGs/TME/", i,"_tumor_vs_tissue_ctrl.csv") )
}

### mouse data tumor vs. colon control 
for(i in cell_types) {
  Idents(mouse) <- "annotation"
  sub <- subset(mouse, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"tumor_wt","adult_colon_wt",FALSE,0.25,
                      paste0("/scratch/khandl/eos_tumor/data_files/DEGs/TME/", i,"_tumor_wt_vs_adult_colon_wt.csv") )
}

### mouse data tumor PHIL vs. WT 
for(i in cell_types) {
  Idents(mouse) <- "annotation"
  sub <- subset(mouse, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"tumor_phil","tumor_wt",FALSE,0.25,
                      paste0("/scratch/khandl/eos_tumor/data_files/DEGs/TME/", i,"_tumor_phil_vs_tumor_wt.csv") )
}

##### plotting of DEGs of interest
### human data tumor vs. tissue ctrl 
## load DEG csv files 
mono <- read.csv("/scratch/khandl/eos_human/data_files/DEGs/TME/Monocytes_tumor_vs_tissue_ctrl.csv")
mac <- read.csv("/scratch/khandl/eos_human/data_files/DEGs/TME/Macrophages_tumor_vs_tissue_ctrl.csv")
tams <- read.csv("/scratch/khandl/eos_human/data_files/DEGs/TME/TAMs_tumor_vs_tissue_ctrl.csv")

## extract genes of interest and add 0 if needed 
genes_of_interest <- c("IL1RN","ADAM19","CCL3","NFKBIA","PLA2G7","ADAM8","IL7R","IL1R1","VEGFA","CXCL2","CD300A","CXCL8","ADAM9",
                       "SELL","S100A8","CD55","HLA-B","HLA-E","JUND","CXCL1","IL6","CD80","CXCL3","TNF","HLA-C","APOE",
                       "LYVE1","CXCL12","F13A1","FOLR2","THBS1","FOSB","SPP1"
)

## prepare df with genes of interest
df_mono <- add_missing_genes_to_df(mono, genes_of_interest,"mono" )
df_mac <- add_missing_genes_to_df(mac, genes_of_interest,"mac" )
df_tams <- add_missing_genes_to_df(tams, genes_of_interest,"tams" )

## merge dataframe 
merged <- merge(df_mono,df_mac, by = "row.names")
rownames(merged) <- merged$Row.names
merged <- merge(merged,df_tams, by = "row.names")
rownames(merged) <- merged$Row.names

merged$Row.names <- NULL
merged$Row.names <- NULL

merged$mono <- as.numeric(merged$mono)
merged$mac <- as.numeric(merged$mac)
merged$tams <- as.numeric(merged$tams)

# plot 
ComplexHeatmap::Heatmap(merged, name=paste0("Scaled log2FC"),
                        column_names_gp = grid::gpar(fontsize = 10),
                        row_names_gp = grid::gpar(fontsize = 10),
                        column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                        column_names_rot = 45,cluster_rows = TRUE, cluster_columns = TRUE,column_title = "DEGs")

### mouse data tumor vs. colon control 
## load DEG csv files 
mono <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Monocytes_tumor_wt_vs_adult_colon_wt.csv")
mac <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Macrophages_tumor_wt_vs_adult_colon_wt.csv")
tams <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/TAMs_tumor_wt_vs_adult_colon_wt.csv")

## genes of interest 
genes_of_interest <- c("Cd36","Cd55","Adam19","F13a1","Alox5","Relb","Nfkb1","H2-K2","Jun","Tgfbr2","Fn1","Arg1","Apoe",
                       "Ccl6","Ccl9","C1qb","Mmp19","Spp1","S100a1","S100a10","Vegfa","Cd72",
                       "Lyz2","Thbs1","Fcgr2b","Tgfbi","H2-K1","Spp1","Fn1","Nfkbia","Nfkbid","H2-K1","H2-D1",
                       "Pla2g7","Cxcr4","Rela","Mki67","Adgre4","Rel","Cd300a",
                       "Il6ra","Il1a","Ly6e","Il1rn","Ccr2","Ccl3","Ccl4","S100a11",
                       "Csf1r","Adora2b","Mertk","Mrc1"
)

## prepare df with genes of interest
df_mono <- add_missing_genes_to_df(mono, genes_of_interest,"mono" )
df_mac <- add_missing_genes_to_df(mac, genes_of_interest,"mac" )
df_tams <- add_missing_genes_to_df(tams, genes_of_interest,"tams" )

## merge dataframe 
merged <- merge(df_mono,df_mac, by = "row.names")
rownames(merged) <- merged$Row.names
merged <- merge(merged,df_tams, by = "row.names")
rownames(merged) <- merged$Row.names

merged$Row.names <- NULL
merged$Row.names <- NULL

merged$mono <- as.numeric(merged$mono)
merged$mac <- as.numeric(merged$mac)
merged$tams <- as.numeric(merged$tams)

# plot 
ComplexHeatmap::Heatmap(merged, name=paste0("Scaled log2FC"),
                        column_names_gp = grid::gpar(fontsize = 10),
                        row_names_gp = grid::gpar(fontsize = 10),
                        column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                        column_names_rot = 45,cluster_rows = TRUE, cluster_columns = TRUE,column_title = "DEGs")

### mouse tumor PHIL vs. WT 
## load DEG csv files 
mono <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Monocytes_tumor_phil_vs_tumor_wt.csv")
mac <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Macrophages_tumor_phil_vs_tumor_wt.csv")
tams <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/TAMs_tumor_phil_vs_tumor_wt.csv")

## genes of interest 
genes_of_interest <- c("Cd36","Cd55","Adam19","F13a1","Alox5","Relb","Nfkb1","H2-K2","Jun","Tgfbr2","Fn1","Arg1","Apoe",
                       "Ccl6","Ccl9","C1qb","Mmp19","Spp1","S100a1","S100a10","Vegfa","Cd72",
                       "Lyz2","Thbs1","Fcgr2b","Tgfbi","H2-K1","Spp1","Fn1","Nfkbia","Nfkbid","H2-K1","H2-D1",
                       "Pla2g7","Cxcr4","Rela","Mki67","Adgre4","Rel","Cd300a",
                       "Il6ra","Il1a","Ly6e","Il1rn","Ccr2","Ccl3","Ccl4","S100a11",
                       "Csf1r","Adora2b","Mertk","Mrc1"
)

## prepare df with genes of interest
df_mono <- add_missing_genes_to_df(mono, genes_of_interest,"mono" )
df_mac <- add_missing_genes_to_df(mac, genes_of_interest,"mac" )
df_tams <- add_missing_genes_to_df(tams, genes_of_interest,"tams" )

## merge dataframe 
merged <- merge(df_mono,df_mac, by = "row.names")
rownames(merged) <- merged$Row.names
merged <- merge(merged,df_tams, by = "row.names")
rownames(merged) <- merged$Row.names

merged$Row.names <- NULL
merged$Row.names <- NULL

merged$mono <- as.numeric(merged$mono)
merged$mac <- as.numeric(merged$mac)
merged$tams <- as.numeric(merged$tams)

# plot 
ComplexHeatmap::Heatmap(merged, name=paste0("Scaled log2FC"),
                        column_names_gp = grid::gpar(fontsize = 10),
                        row_names_gp = grid::gpar(fontsize = 10),
                        column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                        column_names_rot = 45,cluster_rows = TRUE, cluster_columns = TRUE,column_title = "DEGs")

### plot DEGS of TAMs in a barplot 
genes_of_interest <- c("Slpi","Ccl5","Ccl4","Cxcl10","Ccl3","Cxcl1","Cxcl3","Cxcl2","Ptgs2","Klf6","Nlrp3","Il1b","Tnf","C1qb","Tgfbi",
         "Chil3","Mrc1","Spp1","Ccl6","Mertk","Csf2rb2","Arg1","Fn1","Mmp19","Vegfa")

tams <- tams[tams$X %in% genes_of_interest,]

p <- ggplot(tams, aes(x = reorder(X, -avg_log2FC), y = avg_log2FC, fill = p_val_adj)) + geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/eos_tumor/figures/PHIL_wt.svg", width = 12, height = 6, plot = p)

