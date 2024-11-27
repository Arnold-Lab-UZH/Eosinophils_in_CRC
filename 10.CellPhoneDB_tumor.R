########## This code applies CellPhoneDB to predict ligand-receptor interactions in the tumor microenvironment (human and mouse) ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
mouse <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")
human <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")

# extract only tumor samples 
Idents(mouse) <- "condition"
mouse <- subset(mouse, idents = "tumor_wt")
Idents(human) <- "tissue"
human <- subset(human, idents = "tumor")

##### generate input for CellPhoneDB 
Input_files_CellPhoneDB_generation_hs(human, "annotation","tumor", "/scratch/khandl/eos_human/CellPhoneDB/input/") 
Input_files_CellPhoneDB_generation_mm(mouse, "annotation","tumor", "/scratch/khandl/eos_tumor/CellPhoneDB/input/") 

##### run in terminal within a conda environment and python 3
### human 
conda activate cpdb
python3
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = '/mnt/cpdb/v5.0.0/cellphonedb.zip',
  meta_file_path = '/mnt/cpdb/input/tumor_meta.txt',
  counts_file_path = '/mnt/cpdb/input/tumor_count.txt',
  counts_data = 'hgnc_symbol',
  score_interactions = True,
  threshold = 0.1,
  output_path = '/mnt/cpdb/output/human_tumor')

### mouse  
cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = '/mnt/cpdb/v5.0.0/cellphonedb.zip',
  meta_file_path = '/mnt/cpdb/input/tumor_wt_meta.txt',
  counts_file_path = '/mnt/cpdb/input/tumor_wt_count.txt',
  counts_data = 'ensembl',
  score_interactions = True,
  threshold = 0.1,
  output_path = '/mnt/cpdb/output/tumor_wt')

##### plotting of the number of interactions between eos and any other cell type 
### human 
means <- read.delim("/data/khandl/CellPhoneDB/eos_human/human_tumor/statistical_analysis_significant_means_08_13_2024_120238.txt",check.names = FALSE)
means[is.na(means)] <- 0

eos_tams <- ligand_receptor_counts(means, "Eosinophils","TAMs")
eos_cd4 <- ligand_receptor_counts(means, "Eosinophils","CD4_T")
eos_cd8 <- ligand_receptor_counts(means, "Eosinophils","CD8_T")
eos_mast <- ligand_receptor_counts(means, "Eosinophils","Mast")
eos_dc <- ligand_receptor_counts(means, "Eosinophils","DCs")
eos_mono <- ligand_receptor_counts(means, "Eosinophils","Monocytes")
eos_macro <- ligand_receptor_counts(means, "Eosinophils","Macrophages")
eos_pcs <- ligand_receptor_counts(means, "Eosinophils","PC")
eos_matureB <- ligand_receptor_counts(means, "Eosinophils","B_mature")
eos_neut <- ligand_receptor_counts(means, "Eosinophils","Neutrophils")
eos_eos <- ligand_receptor_counts(means, "Eosinophils","Eosinophils")
eos_baso <- ligand_receptor_counts(means, "Eosinophils","Basophils")
eos_fib <- ligand_receptor_counts(means, "Eosinophils","Fibroblasts")
eos_epi <- ligand_receptor_counts(means, "Eosinophils","Epithelial")
eos_endo <- ligand_receptor_counts(means, "Eosinophils","Endothelial")

a <- c("Eos|TAMs","Eos|CD4_T","Eos|CD8_T","Eos|Mast","Eos|DCs","Eos|Mono","Eos|Mac","Eos|PCs","Eos|matureB","Eos|Neutro","Eos|Eos","Eos|Basophils")
b <- c(eos_tams, eos_cd4, eos_cd8,eos_mast,eos_dc, eos_mono,eos_macro,eos_pcs,eos_matureB,eos_neut,eos_eos,eos_baso)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","total_counts")
df$total_counts <- as.numeric(df$total_counts)

p <- ggplot(df, aes(x = interacting_pair, y = total_counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  
ggsave("/scratch/khandl/eos_tumor/figures/LR/LRcounts_human.svg", width = 12, height = 6, plot = p)

### mouse 
means <- read.delim("/data/khandl/CellPhoneDB/eos_tumor/tumor_wt/statistical_analysis_significant_means_05_21_2024_131315.txt",check.names = FALSE)
means[is.na(means)] <- 0

eos_tams <- ligand_receptor_counts(means, "Eosinophils","TAMs")
eos_cd4 <- ligand_receptor_counts(means, "Eosinophils","CD4_T")
eos_cd8 <- ligand_receptor_counts(means, "Eosinophils","CD8_T")
eos_ilc <- ligand_receptor_counts(means, "Eosinophils","ILCs")
eos_dc <- ligand_receptor_counts(means, "Eosinophils","DCs")
eos_mono <- ligand_receptor_counts(means, "Eosinophils","Monocytes")
eos_macro <- ligand_receptor_counts(means, "Eosinophils","Macrophages")
eos_pcs <- ligand_receptor_counts(means, "Eosinophils","PCs")
eos_matureB <- ligand_receptor_counts(means, "Eosinophils","mature_B")
eos_neut <- ligand_receptor_counts(means, "Eosinophils","Neutrophils")
eos_eos <- ligand_receptor_counts(means, "Eosinophils","Eosinophils")
eos_trm <- ligand_receptor_counts(means, "Eosinophils","TRM")

a <- c("Eos|TAMs","Eos|CD4_T","Eos|CD8_T","Eos|ILCs","Eos|DCs","Eos|Mono","Eos|Mac","Eos|PCs","Eos|matureB","Eos|Neutro","Eos|Eos","Eos|TRM")
b <- c(eos_tams, eos_cd4, eos_cd8,eos_ilc,eos_dc, eos_mono,eos_macro,eos_pcs,eos_matureB,eos_neut,eos_eos,eos_trm)

df <- as.data.frame(a,b)
df$counts <- rownames(df)
colnames(df) <- c("interacting_pair","total_counts")
df$total_counts <- as.numeric(df$total_counts)

p <- ggplot(df, aes(x = interacting_pair, y = total_counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))  
ggsave("/scratch/khandl/eos_tumor/figures/LR/LRcounts.svg", width = 12, height = 6, plot = p)

##### comparison of sig. interactions between Eos and SPP1 TAMs in mouse and human tumors 
### mouse
df <- read.delim("/data/khandl/CellPhoneDB/eos_tumor/tumor_wt/statistical_analysis_significant_means_05_21_2024_131315.txt",check.names = FALSE)
df[is.na(df)] <- 0

df <- df[,colnames(df) %in% c("interacting_pair","directionality","Eosinophils|TAMs")]
df2 <- df[which(df$`Eosinophils|TAMs` > 0),]
rownames(df2) <- df2$interacting_pair
df2$directionality <- NULL

df2$celltypes <- "Eosinophils_TAMs"
colnames(df2) <- c("interacting_pair","mean","celltypes")
df_mouse <- df2
mouse_interactions <- df_mouse$interacting_pair

### human
df <- read.delim("/data/khandl/CellPhoneDB/eos_human/human_tumor/statistical_analysis_significant_means_08_13_2024_120238.txt",check.names = FALSE)
df[is.na(df)] <- 0

df <- df[,colnames(df) %in% c("interacting_pair","directionality","Eosinophils|TAMs")]
df2 <- df[which(df$`Eosinophils|TAMs` > 0),]
rownames(df2) <- df2$interacting_pair
df2$directionality <- NULL

df2$celltypes <- "Eosinophils_TAMs"
colnames(df2) <- c("interacting_pair","mean","celltypes")
human_interactions <- df2$interacting_pair

### find the shared interactions 
shared_interactions <- intersect(mouse_interactions, human_interactions)

### plot pvalues and means of shared interactions  
# mouse
pval <- read.delim("/data/khandl/CellPhoneDB/eos_tumor/tumor_wt/statistical_analysis_pvalues_05_21_2024_131315.txt",check.names = FALSE)
pval <- pval[pval$interacting_pair %in% shared_interactions ,colnames(pval) %in% c("interacting_pair","Eosinophils|TAMs")]
means <- read.delim("/data/khandl/CellPhoneDB/eos_tumor/tumor_wt/statistical_analysis_means_05_21_2024_131315.txt",check.names = FALSE)
means <- means[means$interacting_pair %in% shared_interactions ,colnames(means) %in% c("interacting_pair","Eosinophils|TAMs")]

colnames(pval) <- c("interacting_pair","pval")
colnames(means) <- c("interacting_pair","means")
df_mouse <- merge(pval,means, by = "interacting_pair")
df_mouse$celltype <- "TAMs_mouse"

# human
pval <- read.delim("/data/khandl/CellPhoneDB/eos_human/human_tumor/statistical_analysis_pvalues_08_13_2024_120238.txt",check.names = FALSE)
pval <- pval[pval$interacting_pair %in% shared_interactions ,colnames(pval) %in% c("interacting_pair","Eosinophils|TAMs")]
means <- read.delim("/data/khandl/CellPhoneDB/eos_human/human_tumor/statistical_analysis_means_08_13_2024_120238.txt",check.names = FALSE)
means <- means[means$interacting_pair %in% shared_interactions ,colnames(means) %in% c("interacting_pair","Eosinophils|TAMs")]

colnames(pval) <- c("interacting_pair","pval")
colnames(means) <- c("interacting_pair","means")
df_human <- merge(pval,means, by = "interacting_pair")
df_human$celltype <- "TAMs_human"

# merge and plot 
merged <- rbind(df_mouse,df_human)

p <- ggplot(merged, aes(x = means, y = interacting_pair,size = -pval, color = celltype)) + 
  geom_point() + scale_size(name = "padj", range = c(10, 20)) + theme_light() + 
  scale_color_manual(values = c(TAMs_human="#371D3D", TAMs_mouse="#929992"))
ggsave("/scratch/khandl/eos_tumor/figures/lr/eos_tams_human_mouse.svg", width = 12, height = 6, plot = p)




