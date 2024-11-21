########## This code runs SCENIC analysis of human eosinophil clusters after integration ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
obj <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated.rds")

##### download the motive databases from here  
#https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/
#hg19-500bp-upstream-7species.mc9nr.feather
#hg19-tss-centered-10kb-7species.mc9nr.feather 

setwd("/data/khandl/SCENIC/human_eos/SCENIC_clusters") 

##### downsample to 500 per cluster to speed up analysis 
Idents(obj) <- "mnn.clusters"
set.seed(111)
sub <- subset(x = obj, downsample = 500)

##### run SCENIC analysis 
### Initialize settings
exprMat <- as.matrix(sub[["RNA"]]$counts)
cellInfo <- sub@meta.data

colnames(cellInfo)[which(colnames(cellInfo)=="mnn.clusters")] <- "CellType"
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

org="hgnc" 
dbDir="/home/khandl/data/common/SCENIC/human" # RcisTarget databases location
myDatasetTitle="SCENIC human eos clusters" 
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs

### load motive annotation and rename it 
data(list="motifAnnotations_hgnc", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### filter genes 
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

### calculate correlation
runCorrelation(exprMat_filtered, scenicOptions)

### run GENIE3: infer potential transcription factor targets based on the expression data
exprMat_filtered <- log2(exprMat_filtered+1) 
motifAnnotations_mgi_v8 <- motifAnnotations
runGenie3(exprMat_filtered, scenicOptions)

### score GRN (regulons) 
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

exprMat_log <- log2(exprMat+1)
dim(exprMat)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipHeatmap = TRUE, skipTsne = TRUE)

### Binarizing the network
#to make it on/off, makes it more clear to compare between conditions 
#Clustering / dim reduction on the regulon activity.
nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them

runSCENIC_4_aucell_binarize(scenicOptions, skipBoxplot = TRUE, skipHeatmaps = TRUE, skipTsne = TRUE)

### save output 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### plot regulons in heatmap
scenicOptions <- readRDS(file="int/scenicOptions.Rds") 

### extract activity scores per regulon and cluster and scale 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

### save regulon outputs in csv files 
write.csv(regulonActivity_byCellType_Scaled, "/data/khandl/SCENIC/human_eos/SCENIC_clusters/Human_eos_clusters_regulons_scaled_actvity.csv")
write.csv(regulonActivity_byCellType, "/data/khandl/SCENIC/human_eos/SCENIC_clusters/Human_eos_clusters_regulons_actvity.csv")

### plot regulons of interest in heatmap 
regulonActivity_byCellType <- read.csv("/data/khandl/SCENIC/human_eos/SCENIC_clusters/Human_eos_clusters_regulons_actvity.csv")
rownames(regulonActivity_byCellType) <- regulonActivity_byCellType$X
regulonActivity_byCellType$X <- NULL

# select regulons of interest 
regulons_specific <- c("FOSB (23g)", "FOSL2 (59g)","NFKB2 (57g)","RELB (172g)","NFKB1 (97g)","REL (62g)" ,
                       "RARA_extended (38g)","IRF1_extended (61g)","FOSL1_extended (187g)"

)

regulonActivity_byCellType2 <- regulonActivity_byCellType[rownames( regulonActivity_byCellType) %in% regulons_specific,]
regulonActivity_byCellType2$X8 <- NULL
regulonActivity_byCellType2$X9 <- NULL

# plot 
ComplexHeatmap::Heatmap(regulonActivity_byCellType2, name=paste0("regulon activity"),
                        column_names_gp = grid::gpar(fontsize = 10),
                        row_names_gp = grid::gpar(fontsize = 10),
                        column_dend_height = unit(0.6, "cm"), row_dend_width = unit(0.85, "cm"),
                        column_names_rot = 45,cluster_rows = TRUE, cluster_columns = TRUE,column_title = "Regulons")


