########## This code applies CytoSig analysis to mouse eosinophils ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
human_data <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")
mouse_data <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")
human_data <- JoinLayers(human)
mouse_data <- JoinLayers(mouse)

### load human and mouse ensemble symbols
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

##### mouse data 
### Input generation between tumor and control 
input_cytosig_mouse_2cond(mouse_data, "condition","tumor_wt","adult_colon_wt",0.1, "tumor_wt_vs_adult_colon_wt",
                          "/scratch/khandl/eos_tumor/chemo_cyto/input_cytosig_tumor_wt_vs_adult_colon_wt.csv") 

### run on the web-based tool https://cytosig.ccr.cancer.gov/
## upload outputs from https://cytosig.ccr.cancer.gov/profiler/task_run/ to the server 
#scp /Users/handler/Downloads/cytosig_tumor_wt_vs_adult_colon_wt.csv khandl@cluster.s3it.uzh.ch:data/eos_human/cytosig/ 

### plot significant CytoSig results 
file <- read.csv("/data/khandl/eos_human/cytosig/cytosig_tumor_wt_vs_adult_colon_wt.csv")
file <- file[,c(1,4,5)]

#remove IDs that are non-significant in both conditions 
file <- file[file$Pvalue <= 0.05,]
print(file)

rownames(file) <- file$ID
file$ID <- NULL
file$Pvalue <- NULL

print(ComplexHeatmap::Heatmap(file, name="tumor_vs_tissue_ctrl",
                              column_names_gp = grid::gpar(fontsize = 10),
                              row_names_gp = grid::gpar(fontsize = 5),
                              column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                              column_names_rot = 45,cluster_rows = TRUE, cluster_columns = FALSE,column_title = "GO terms"))

##### human data 
### Input generation between tumor and control 
input_cytosig_human_2cond(human_data, "tissue","tumor","tissue_ctrl",0.1, "tumor_vs_tissue_ctrl",
                          "/scratch/khandl/eos_human/figures/chemo_cyto/input_cytosig_tumor_vs_tissue_ctrl.csv") 

### run on the web-based tool https://cytosig.ccr.cancer.gov/
## upload outputs from https://cytosig.ccr.cancer.gov/profiler/task_run/ to the server 
#scp /Users/handler/Downloads/cytosig_results_tumor_vs_tissue_ctrl.csv khandl@cluster.s3it.uzh.ch:data/eos_human/cytosig/ 

### plot significant CytoSig results 
file <- read.csv("/data/khandl/eos_human/cytosig/cytosig_results_tumor_vs_tissue_ctrl.csv")
file <- file[,c(1,4,5)]

#remove IDs that are non-significant in both conditions 
file_isg <- file[file$Pvalue <= 0.05,]
print(file)

file <- file[file$ID %in% c(file_isg$ID,"IL33","IL10","HMGB1","TGFB3","IL6","IL18","OSM","IL1B","IL1A","BMP4","CXCL12","IL15","CCL2","WNT3A","TGFB1"),]

rownames(file) <- file$ID
file$ID <- NULL
file$Pvalue <- NULL

print(ComplexHeatmap::Heatmap(file, name="tumor_vs_tissue_ctrl",
                              column_names_gp = grid::gpar(fontsize = 10),
                              row_names_gp = grid::gpar(fontsize = 5),
                              column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                              column_names_rot = 45,cluster_rows = TRUE, cluster_columns = FALSE,column_title = "GO terms"))

