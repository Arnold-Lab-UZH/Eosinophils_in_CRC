########## This code compares eosinophils between human and mouse ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
human <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/patients_colonic_annotated.rds")
mouse <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_colonic_annotated.rds")
human <- JoinLayers(human)
mouse <- JoinLayers(mouse)

##### human data 
Idents(human) <- "annotation"
p <- DotPlot(human, features = c("TGFB1","IL10"))
ggsave("/scratch/khandl/eos_human/tgfb1_il10/human_data_dotplot.svg", width = 12, height = 6, plot = p)

Idents(human) <- "annotation"
p <- DotPlot(human, features = c("IL33","HMGB1"))
ggsave("/scratch/khandl/eos_human/tgfb1_il10/human_data_dotplot_IL33_HMGB1.svg", width = 12, height = 6, plot = p)

## tumor vs. NAT
for (i in cell_types) {
  Idents(human) <- "annotation"
  sub <- subset(human, idents = i)
  Idents(sub) <- "tissue"
  DEG_to_csv_two_cond(sub,"tumor","tissue_ctrl",FALSE,0.25,
                      paste0("/scratch/khandl/eos_tumor/data_files/DEGs/TME/", i,"_tumor_vs_NAT.csv") )
}

#extract Tgfb1 and Il10 
monocytes <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Monocytes_tumor_vs_NAT.csv")
macrophages <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Macrophages_tumor_vs_NAT.csv")
TAMs <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/TAMs_tumor_vs_NAT.csv")

#extract Tgfb1 and Il10 
monocytes <- monocytes[monocytes$X %in% c("Tgfb1","Il10"),]
macrophages <- macrophages[macrophages$X %in% c("Tgfb1","Il10"),]
TAMs <- TAMs[TAMs$X %in% c("Tgfb1","Il10"),]

##### mouse data 
Idents(mouse) <- "annotation"
p <- DotPlot(mouse, features = c("Tgfb1","Il10"))
ggsave("/scratch/khandl/eos_human/tgfb1_il10/mouse_data_dotplot.svg", width = 12, height = 6, plot = p)

Idents(mouse) <- "annotation"
p <- DotPlot(mouse, features = c("Il33","Hmgb1"))
ggsave("/scratch/khandl/eos_human/tgfb1_il10/mouse_data_dotplot_il33_hmgb1.svg", width = 12, height = 6, plot = p)

### DEG analysis 
cell_types <- c("TAMs","Monocytes","Macrophages")

## tumor vs. control 
for (i in cell_types) {
  Idents(mouse) <- "annotation"
  sub <- subset(mouse, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"tumor_wt","adult_colon_wt",FALSE,0.25,
                      paste0("/scratch/khandl/eos_tumor/data_files/DEGs/TME/", i,"_tumor_wt_vs_control.csv") )
}

## tumor phil vs. wt 
for (i in cell_types) {
  Idents(mouse) <- "annotation"
  sub <- subset(mouse, idents = i)
  Idents(sub) <- "condition"
  DEG_to_csv_two_cond(sub,"tumor_phil","tumor_wt",FALSE,0.25,
                      paste0("/scratch/khandl/eos_tumor/data_files/DEGs/TME/", i,"_tumor_phil_vs_tumor_wt.csv") )
}

#extract Tgfb1 and Il10 
monocytes <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Monocytes_tumor_phil_vs_tumor_wt.csv")
macrophages <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/Macrophages_tumor_phil_vs_tumor_wt.csv")
TAMs <- read.csv("/scratch/khandl/eos_tumor/data_files/DEGs/TME/TAMs_tumor_phil_vs_tumor_wt.csv")

#extract Tgfb1 and Il10 
monocytes <- monocytes[monocytes$X %in% c("Tgfb1","Il10"),]
macrophages <- macrophages[macrophages$X %in% c("Tgfb1","Il10"),]
TAMs <- TAMs[TAMs$X %in% c("Tgfb1","Il10"),]
