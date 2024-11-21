########## This code analyses DEGs of Spp1 high TAMS between tumor PHIL vs. WT  ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
mouse <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")

##### run DEG analysis 
Idents(mouse) <- "annotation"
sub <- subset(mouse, idents = "TAMs")
Idents(sub) <- "condition"
DEG_to_csv_two_cond(sub,"tumor_wt","adult_colon_wt",FALSE,0.25,
                    paste0("/scratch/khandl/eos_tumor/data_files/DEGs/TME/", "TAMs","_tumor_phil_vs_tumor_wt.csv") )


##### plotting of DEGs of interest
### load DEG csv files 
tams <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/TAMs_tumor_phil_vs_tumor_wt.csv")

### plot DEGS of TAMs in a barplot 
genes_of_interest <- c("Slpi","Ccl5","Ccl4","Cxcl10","Ccl3","Cxcl1","Cxcl3","Cxcl2","Ptgs2","Klf6","Nlrp3","Il1b","Tnf","C1qb","Tgfbi",
         "Chil3","Mrc1","Spp1","Ccl6","Mertk","Csf2rb2","Arg1","Fn1","Mmp19","Vegfa")

tams <- tams[tams$X %in% genes_of_interest,]

p <- ggplot(tams, aes(x = reorder(X, -avg_log2FC), y = avg_log2FC, fill = p_val_adj)) + geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90)) 
ggsave("/scratch/khandl/eos_tumor/figures/PHIL_wt.svg", width = 12, height = 6, plot = p)

