########## This code analyses the MHC-I and antigen processing score in eosinophils between tumor and control (human and mouse) ##########
### use the DEG analysis results from 7.Eosinophils_DEG_analysis.R

##### link to libraries and functions
source("~/Projects/Eosinophils_in_CRC/1.Packages_and_functions.R")

##### load R object 
mouse_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_eosinophils_integrated_annotated.rds")
human_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated_annotated.rds")

##### mouse data 
MHCI_antigen_processing_vln_2_cond_mm(mouse_eos, "condition",c("adult_colon_wt", "tumor_wt"),
                                      "condition","adult_colon_wt","tumor_wt", 
                                      c("#848BA3","#29489C"),c(-0.5,3),
                                      "/scratch/khandl/eos_tumor/figures/signature_scores/MHCI_adult_colon_wt_vs_tumor_wt.svg")

##### human data 
### tumor vs. tissue ctrl 
MHCI_antigen_processing_vln_2_cond_hs(human_eos, "tissue",c("tissue_ctrl", "tumor"),
                                      "tissue","tissue_ctrl","tumor", 
                                      c("#9E5E6A","#8A181A"),c(-0.5,3),
                                      "/scratch/khandl/eos_human/figures/signature_scores/MHCI_control_vs_tumor.svg")

### cluster 2 vs. cluster 1
MHCI_antigen_processing_vln_2_cond_hs(human_eos, "mnn.clusters",c("1", "2"),
                                      "mnn.clusters","1","2", 
                                      c("#CA13EA","#EAA10F"),c(-0.5,3),
                                      "/scratch/khandl/eos_human/figures/signature_scores/MHCI_cluster1_vs_cluster2.svg")
