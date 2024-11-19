##########  Pre-processing of human scRNAseq data ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

#use the standard first cutoff settings, min.cell =3, min.features = 200

### P1 
P1_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P1_tumor_Expression_Data.st"), 
                                                      "P1",3,200,  "P1_tumor","tumor","Exp1","patient")
P1_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "human_fixed_counts", "P1_NAT_Expression_Data.st"), 
                                                        "P1",3,200, "P1_tissue_ctrl","tissue_ctrl","Exp1","patient")

### P2 
P2_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P2_tumor_Expression_Data.st"), 
                                                      "P2",3,200,  "P2_tumor","tumor","Exp2","patient")
P2_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P2_NAT_Expression_Data.st"), 
                                                        "P2",3,200,  "P2_tissue_ctrl","tissue_ctrl","Exp2","patient")

### P3
P3_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P3_tumor_Expression_Data.st"), 
                                                      "P3",3,200,  "P3_tumor","tumor","Exp3","patient")
P3_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P3_NAT_Expression_Data.st"), 
                                                        "P3",3,200,  "P3_tissue_ctrl","tissue_ctrl","Exp3","patient")

### P4
P4_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P4_tumor_Expression_Data.st"), 
                                                      "P4",3,200,  "P4_tumor","tumor","Exp4","patient")
P4_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P4_NAT_Expression_Data.st"), 
                                                        "P4",3,200,  "P4_tissue_ctrl","tissue_ctrl","Exp4","patient")

### P5
P5_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P5_tumor_Expression_Data.st"), 
                                                      "P5",3,200,  "P5_tumor","tumor","Exp5","patient")
P5_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P5_NAT_Expression_Data.st"), 
                                                        "P5",3,200,  "P5_tissue_ctrl","tissue_ctrl","Exp5","patient")

### P6
P6_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P6_tumor_Expression_Data.st"), 
                                                      "P6",3,200,  "P6_tumor","tumor","Exp7","patient")
P6_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P6_NAT_Expression_Data.st"), 
                                                        "P6",3,200,  "P6_tissue_ctrl","tissue_ctrl","Exp7","patient")

### P7 
P7_tumor <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P7_tumor_Expression_Data.st"), 
                                                      "P7",3,200,  "P7_tumor","tumor","Exp8","patient")
P7_control <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "P7_NAT_Expression_Data.st"), 
                                                        "P7",3,200,  "P7_tissue_ctrl","tissue_ctrl","Exp8","patient")
### healthy individuals 
H1 <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "Blood_healthy1_hs_Expression_Data.st"), 
                                                "H1",3,200,  "H1_blood","blood_healthy","Exp6","healthy")
H2 <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "Blood_healthy2_hs_Expression_Data.st"), 
                                                "H2",3,200,  "H2_blood","blood_healthy","Exp6","healthy")
H3 <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "Blood_healthy3_hs_Expression_Data.st"), 
                                                "H3",3,200,  "H3_blood","blood_healthy","Exp6","healthy")
H4 <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "Blood_healthy4_hs_Expression_Data.st"), 
                                                "H4",3,200,  "H4_blood","blood_healthy","Exp6","healthy")
H5 <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "Blood_healthy5_hs_Expression_Data.st"), 
                                                "H5",3,200,  "H5_blood","blood_healthy","Exp6","healthy")
H6 <- create_seurat_plus_DecontX_human_biopsies(file.path("/data/khandl/raw_data/", "counts", "Blood_healthy6_hs_Expression_Data.st"), 
                                                "H6",3,200,  "H6_blood","blood_healthy","Exp6","healthy")

### Merge samples
patients <- merge(P1_tumor, y = c(P1_control, P2_tumor,P2_control,P3_tumor,P3_control,P4_tumor,P4_control,P5_tumor,P5_control,P6_tumor,P6_control,P7_tumor,P7_control),
                  add.cell.ids = c("h1", "h2","h3","h4", "h5","h6","h7", "h8","h9","h10", "h11","h12", "h13", "h14"))
patients <- JoinLayers(patients)

blood_healthy <- merge(H1, y = c(H2,H3,H4,H5,H6),
                       add.cell.ids = c("b1", "b2","b3","b4", "b5","b6"))
blood_healthy <- JoinLayers(blood_healthy)

### Add mitochondrial percentage per cell and apply appropriate quality cutoffs 
patients$percent.mt <- PercentageFeatureSet(patients, pattern = "^MT.")
blood_healthy$percent.mt <- PercentageFeatureSet(blood_healthy, pattern = "^MT.")

### QC plots 
hs <- merge(patients, blood_healthy)

Idents(hs) <- "tissue"
p <- VlnPlot(hs, features = c("nFeature_RNA"), pt.size = 0)
ggsave("/scratch/khandl/eos_human/QC/nFeature.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(hs, features = c("percent.mt"), pt.size = 0)
ggsave("/scratch/khandl/eos_human/QC/mito.svg", width = 8, height = 5, plot = p)
p <- VlnPlot(hs, features = c("nCount_RNA"), pt.size = 0)
ggsave("/scratch/khandl/eos_human/QC/nCount.svg", width = 8, height = 5, plot = p)

### apply mitochondrial and nFeature cutoffs
patients <- subset(patients, subset = percent.mt < 25)
blood_healthy <- subset(blood_healthy, subset = percent.mt < 25)

### save object
saveRDS(patients, file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_combined.rds")
saveRDS(blood_healthy, file = "/data/khandl/Eosinophils_in_CRC/seurat_objects/blood_healthy_combined.rds")


