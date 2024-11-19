##########  Pre-processing of human scRNAseq data ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

#use the standard first cutoff settings, min.cell =3, min.features = 200

### WT 
tumor_wt <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data", "M25", "Tumor_mm_WT_Expression_Data.st"), 
  project = "tumor_wt", condition = "tumor_wt",3,200)

disseminated_wt <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data","M25", "Disseminated_mm_WT_Expression_Data.st"), 
  project = "disseminated_wt", condition = "dis_wt",3,200)

adult_colon_wt <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data","M25", "Colon_mm_WT_Expression_Data.st"), 
  project = "adult_colon_wt", condition = "adult_colon_wt",3,200)

blood_wt <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data","M25", "Blood_mm_no_tumor_WT_Expression_Data.st"), 
  project = "blood_wt", condition = "blood_wt",3,200)

### PHIL 
tumor_phil <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data","M25", "Tumor_mm_PHIL_Expression_Data.st"), 
  project = "tumor_phil", condition = "tumor_phil",3,200)

disseminated_phil <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data","M25", "Disseminated_mm_PHIL_Expression_Data.st"), 
  project = "disseminated_phil", condition = "disseminated_phil",3,200)

adult_colon_phil <- create_seurat_plus_DecontX(
  path_to_st_file = file.path("/data/khandl/raw_data","M25", "Colon_mm_PHIL_Expression_Data.st"), 
  project = "adult_colon_phil", condition = "adult_colon_phil",3,200)

### Merge samples
tumor <- merge(tumor_wt, y = c(disseminated_wt,adult_colon_wt, tumor_phil, disseminated_phil, adult_colon_phil,blood_wt),
               add.cell.ids = c("tumor_wt","disseminated_wt","adult_colon_wt", "tumor_phil", "disseminated_phil", "adult_colon_phil","blood_wt"))
tumor <- JoinLayers(tumor)

### Add mitochondrial percentage per cell 
tumor$percent.mt <- PercentageFeatureSet(tumor, pattern = "^mt.")

### QC plots 
Idents(tumor) <- "condition"
sub <- subset(tumor, idents =c("tumor_wt","adult_colon_wt","blood_wt"))

p <- VlnPlot(sub, features = c("nFeature_RNA"), pt.size = 0)
ggsave("/scratch/khandl/eos_tumor/QC/nFeature.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(sub, features = c("percent.mt"), pt.size = 0)
ggsave("/scratch/khandl/eos_tumor/QC/mito.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(sub, features = c("nCount_RNA"), pt.size = 0)
ggsave("/scratch/khandl/eos_tumor/QC/nCount.svg", width = 8, height = 5, plot = p)

### apply mitochondrial and nFeature cutoffs
tumor <- subset(tumor, subset = nFeature_RNA < 5000 & percent.mt < 25)

### save object
saveRDS(tumor, file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/tumor200_5000_25.rds")

