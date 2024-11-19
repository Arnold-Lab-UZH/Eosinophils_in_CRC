########## This code applies GSEA analysis between tumor and control of eosinophils (human and mouse)  ##########
### use the DEG analysis results from 7.Eosinophils_DEG_analysis.R

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
mouse_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_eosinophils_integrated_annotated.rds")
human_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated_annotated.rds")

##### mouse data 
### load GO terms 
#Biological process 
m_df      <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
BP        <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#Hallmarks
m_df      <- msigdbr(species = "Mus musculus", category = "H")
Hallmarks       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#Reactome
m_df      <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
Reactome       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#KEGG
m_df      <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
KEGG       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### run GSEA 
ranks <- ranks_for_fGSEA("/scratch/khandl/eos_tumor/eosinophils/DEGs/mouse_eos_tumor_wt_vs_control_wt.csv")
fGSEA_to_csv(ranks, "/scratch/khandl/eos_tumor/eosinophils/GSEA/mouse_eos_tumor_wt_vs_control_wt.csv")

### plot GO terms of interest 
df <- read.csv("/scratch/khandl/eos_tumor/eosinophils/GSEA/mouse_eos_tumor_wt_vs_control_wt.csv")
#extract terms oi 
df <- df[df$pathway %in% c("GOBP_CELL_CHEMOTAXIS","GOBP_RESPONSE_TO_BACTERIUM","GOBP_CYTOKINE_PRODUCTION",
                           "GOBP_LEUKOCYTE_CHEMOTAXIS","GOBP_MYELOID_LEUKOCYTE_MIGRATION","HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","GOBP_INFLAMMATORY_RESPONSE","GOBP_DEFENSE_RESPONSE"),]

p <- ggplot(df, aes(x = NES, y = pathway, size = padj,color = padj)) + 
  geom_point() + scale_size(name = "padj", range = c(3, 18)) + theme_light()
ggsave("/scratch/khandl/eos_tumor/figures/GO/eos_colon_tumor_Go.svg", width = 12, height = 6, plot = p)

##### human data 
### load GO terms 
#Biological process 
m_df      <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
BP        <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#Hallmarks
m_df      <- msigdbr(species = "Homo sapiens", category = "H")
Hallmarks       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#Reactome
m_df      <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
Reactome       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
#KEGG
m_df      <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
KEGG       <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### tumor vs. tissue ctrl 
## run GSEA 
ranks <- ranks_for_fGSEA("/scratch/khandl/eos_human/data_files/DEGs/eos_cond/tumor_vs_tissue_ctrl.csv")
fGSEA_to_csv(ranks, "/scratch/khandl/eos_human/data_files/GSEA/human_eos_tumor_wt_vs_control_wt.csv")

## plot GO terms of interest 
df <- read.csv("/scratch/khandl/eos_human/data_files/GSEA/human_eos_tumor_wt_vs_control_wt.csv")

p <- plotEnrichment(KEGG$KEGG_CHEMOKINE_SIGNALING_PATHWAY, ranks)
ggsave(file = "/scratch/khandl/eos_human/figures/chemo_cyto/fGSEA.svg", width = 8, height = 6, plot = p)

### cluster 2 vs. cluster 1
## run GSEA 
ranks <- ranks_for_fGSEA("/scratch/khandl/eos_human/data_files/DEGs/eos_cond/cluster2_vs_cluster1.csv")
fGSEA_to_csv(ranks, "/scratch/khandl/eos_human/data_files/GSEA/human_eos_cluster2_vs_cluster1.csv")

## plot GO terms of interest 
df <- read.csv("/scratch/khandl/eos_human/data_files/GSEA/human_eos_cluster2_vs_cluster1.csv")

#extract terms oi 
df <- df[df$pathway %in% c("GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY","GOBP_INFLAMMATORY_RESPONSE","GOBP_NIK_NF_KAPPAB_SIGNALING",
                           "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),]

p <- ggplot(df, aes(x = NES, y = pathway, size = padj,color = padj)) + 
  geom_point() + scale_size(name = "padj", range = c(3, 18)) + theme_light()
ggsave("/scratch/khandl/eos_human/figures/GO_cluster2_vs_cluster1.svg", width = 12, height = 6, plot = p)


