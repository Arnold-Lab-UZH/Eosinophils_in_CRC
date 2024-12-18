########## This code analyses DEGs betweeen tumor and control of eosinophils (human and mouse)  ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
mouse_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_eosinophils_integrated_annotated.rds")
human_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated_annotated.rds")

##### human data 
### DEG analysis tumor vs. NAT
## run DEG analysis only for Eos_H2
Idents(human_eos) <- "mnn.clusters"
H2_eos <- subset(human_eos, idents = "Eos_H2")
Idents(H2_eos) <- "tissue"
DEG_to_csv_two_cond(H2_eos,"tumor","tissue_ctrl",FALSE,0.25,
                    "/scratch/khandl/eos_human/data_files/DEGs/eos_cond/Eos_H2_tumor_vs_tissue_ctrl.csv") 

## plot DEGs results in volcano plot with written genes of interest 
volcano_DGE_showing_goi("/scratch/khandl/eos_human/data_files/DEGs/eos_cond/Eos_H2_tumor_vs_tissue_ctrl.csv",
                        0.25,0.05,"tumor","tissue_ctrl",c("CCL3L1","CCL3","CCL4L2","CCL4","CD69","HLA-DRA","CXCR4","CXCL8",
                                                          "HLA-C","TNFAIP3","NFKBIA","B2M","ALOX5AP"),
                        c("#7F7C7D","#EB1D43","#EB1D43"),c(-0.25,0.25) ,"/scratch/khandl/eos_human/figures/DEGs/volcano_TUMOR_vs_tissue_ctrl.svg") 

## run DEG analysis for all eos together 
Idents(human_eos) <- "tissue"
DEG_to_csv_two_cond(human_eos,"tumor","tissue_ctrl",FALSE,0.25,
                    "/scratch/khandl/eos_human/data_files/DEGs/eos_cond/tumor_vs_tissue_ctrl.csv") 

## run DEG analysis only for Eos_H1  
Idents(human_eos) <- "mnn.clusters"
H1_eos <- subset(human_eos, idents = "Eos_H1")
Idents(H1_eos) <- "tissue"
DEG_to_csv_two_cond(H1_eos,"tumor","tissue_ctrl",FALSE,0.25,
                    "/scratch/khandl/eos_human/data_files/DEGs/eos_cond/Eos_H1_tumor_vs_tissue_ctrl.csv") 

### DEG analysis Eos_H2 vs. Eos_H1
## run DEG analysis 
Idents(human_eos) <- "mnn.clusters"
DEG_to_csv_two_cond(human_eos,"2","1",FALSE,0.25,
                    "/scratch/khandl/eos_human/data_files/DEGs/eos_cond/Eos_H2_vs_Eos_H1.csv") 

### Venn diagram comparing DEGs up in tumor from all eos, Eos_H1 and Eos_H2
Eos_H2 <- read.csv("/scratch/khandl/eos_human/data_files/DEGs/eos_cond/Eos_H2_tumor_vs_tissue_ctrl.csv")
all_eos <- read.csv("/scratch/khandl/eos_human/data_files/DEGs/eos_cond/tumor_vs_tissue_ctrl.csv")
Eos_H1 <- read.csv("/scratch/khandl/eos_human/data_files/DEGs/eos_cond/Eos_H1_tumor_vs_tissue_ctrl.csv")

#extract significant ones 
Eos_H2 <- Eos_H2[Eos_H2$p_val_adj <= 0.05 &Eos_H2$avg_log2FC >0,]
all_eos <- all_eos[all_eos$p_val_adj <= 0.05 & all_eos$avg_log2FC >0,]
Eos_H1 <- Eos_H1[Eos_H1$p_val_adj <= 0.05 & Eos_H1$avg_log2FC >0,]

x <- list("H2" = Eos_H2$X, "all" =all_eos$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

x <- list("all" =all_eos$X, "H1"=Eos_H1$X)
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

##### mouse data 
### run DEG analysis all eos combined tumor vs. control 
mouse_eos <- JoinLayers(mouse_eos)
Idents(mouse_eos) <- "condition"
DEG_to_csv_two_cond(mouse_eos,"tumor_wt","adult_colon_wt",FALSE,0.25,"/scratch/khandl/eos_tumor/eosinophils/DEGs/mouse_eos_tumor_wt_vs_control_wt.csv") 

## plot DEGs results in volcano plot with written genes of interest 
volcano_DGE_showing_goi("/scratch/khandl/eos_tumor/eosinophils/DEGs/mouse_eos_tumor_wt_vs_control_wt.csv",
                        0.25,0.05,"tumor","colon",c("Ccl3","Cxcr4","CXCL8","Camk1d","B2m","Alox5ap","Thbs1","Ccl6","S100a11",
                                                    "H2-K1","Mmp9","H2-D1","Alox15","Il1b","Il1rn","Ly6e","Retnlg","Vegfa","Nkfbid","Il1r2",
                                                    "H2-Q7","Il1a","Cxcl2","Retnla"),
                        c("#EB1D43","#7F7C7D","#EB1D43"),c(-0.25,0.25) ,"/scratch/khandl/eos_tumor/volcano_TUMOR_vs_control.pdf") 

