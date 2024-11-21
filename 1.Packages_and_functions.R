# Packages 
library(Seurat)
library(dplyr)
library(tidyr)
library(decontX)
library(Matrix)
library(ggplot2)
library(SeuratWrappers)
library(pheatmap)
library(ComplexHeatmap)
library(SingleR)
library(celldex)
library(orthogene)
library(scProportionTest)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(biomaRt)
library(ggVennDiagram)
library(RColorBrewer)
library(SCENIC)
library(GENIE3)

##### Functions for pre-processing 
### Function to read in data generated from BD Seven Bridges platform (by Simona Baghai)
data_to_sparse_matrix <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Bioproduct, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Bioproduct")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}

### Function to generate Seurat objects for each sample + application of DecontX (accounting for cell free RNA contamination) 
## for mouse data 
create_seurat_plus_DecontX <- function(
    path_to_st_file,
    project,
    condition, 
    min.cells,
    min.features
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  sce <- as.SingleCellExperiment(condition_sample)
  # run decontX with default settings 
  sce.delta <- decontX(sce)
  # convert back to a Seurat object 
  seuratObject <- CreateSeuratObject(round(decontXcounts(sce.delta)))
  seuratObject$condition <- condition
  
  return(seuratObject)
}

## for human data 
create_seurat_plus_DecontX_human_biopsies <- function(
    path_to_st_file,
    project,
    min.cells,
    min.features,
    condition_oi,
    tissue_oi,
    experiment_oi,
    phenotype_oi
) {
  input_matrix <- data_to_sparse_matrix(path_to_st_file)
  condition_sample <-CreateSeuratObject(input_matrix, 
                                        project = project,
                                        min.cells = min.cells,
                                        min.features = min.features)
  sce <- as.SingleCellExperiment(condition_sample)
  # run decontX with default settings 
  sce.delta <- decontX(sce)
  # convert back to a Seurat object 
  seurat_object <- CreateSeuratObject(round(decontXcounts(sce.delta)))
  seurat_object$condition <- condition_oi
  seurat_object$tissue <- tissue_oi
  seurat_object$experiment <- experiment_oi
  seurat_object$phenotype <- phenotype_oi
  return(seurat_object)
}

##### Functions for cell type annotations 
### Function to project the cell type from SingleR result to the umap space to identify which cluster represents which cell type 
## after FastMNN integration 
project_annotation_to_umap_fastMNN <- function(cell.type, singleResult, seurat_object) {
  temp <- as.data.frame(singleResult[4])
  # singleResult[5] is extracting pruned.labels from data frame
  temp$cell <- rownames(temp)
  temp <- temp%>%filter(temp$pruned.labels %in% cell.type)
  temp <- temp$cell
  print(DimPlot(seurat_object, reduction = "umap.mnn", label = TRUE, label.size = 10, pt.size = 2, cells.highlight = temp, sizes.highlight = 2,raster=FALSE) + NoLegend() + ggtitle(cell.type))
}

## without FastMNN integration 
project_annotation_to_umap <- function(cell.type, singleResult, seurat_object) {
  temp <- as.data.frame(singleResult[4])
  # singleResult[5] is extracting pruned.labels from data frame
  temp$cell <- rownames(temp)
  temp <- temp%>%filter(temp$pruned.labels %in% cell.type)
  temp <- temp$cell
  print(DimPlot(seurat_object, reduction = "umap", label = TRUE, label.size = 10, pt.size = 0.5, cells.highlight = temp, sizes.highlight = 2) + NoLegend() + ggtitle(cell.type))
}

##### Functions for plotting 
### Function to plot genes of interest per condition of interest, scaled per row 
heatmap_goi_coi <- function(
    seurat_object,
    condition_oi,
    markers_oi,
    groups_of_markers,
    number_of_markers_per_group,
    colors_per_group,
    groups_and_colors,
    cluster_rows_cond,
    cluster_cols_cond
){
  # average expression per cluster and condition 
  average_expression <- AverageExpression(seurat_object, return.seurat = FALSE, features = markers_oi, normalization.method = "LogNormalize",assays = "RNA", group.by = condition_oi)
  average_expression_df <- as.data.frame(average_expression)
  average_expression_df <- average_expression_df[match(markers_oi, rownames(average_expression_df)),]
  
  # prepare palette for pheatmap
  paletteLength   <- 50
  myColor         <- colorRampPalette(c("blue", "white", "darkorange"))(paletteLength)
  breaksList      = seq(-2, 2, by = 0.04)
  
  # prepare annotation for pheatmap
  annotation_rows             <- data.frame(markers = rep(groups_of_markers, number_of_markers_per_group))
  rownames(annotation_rows)   <- rownames(average_expression_df)
  annotation_rows$markers     <- factor(annotation_rows$markers, levels = groups_of_markers)
  
  mycolors <- colors_per_group
  names(mycolors) <- unique(annotation_rows$markers)
  mycolors <- list(category = mycolors)
  annot_colors=list(markers=groups_and_colors)
  
  p <- pheatmap(average_expression_df,scale = "row",
                color = colorRampPalette(c("blue", "white", "darkorange"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                breaks = breaksList,
                cluster_rows = cluster_rows_cond, cluster_cols = cluster_cols_cond, 
                border_color = "black", 
                legend_breaks = -2:2, 
                cellwidth = 10, cellheight = 5,
                angle_col = "45", 
                annotation_colors = annot_colors,
                annotation_row = annotation_rows,
                fontsize = 5)
  print(p)
}

##### Functions to compare cell type proportions 
### create table 
create_table_cell_type_prop <- function(
    seurat_object, 
    ident_of_interest1,
    ident_of_interest2,
    output_file_path,
    sample_name
) {
  ann_tab<-  table(seurat_object@meta.data[[ident_of_interest1]], seurat_object@meta.data[[ident_of_interest2]])
  ann_tab <- cbind(ann_tab, Total = rowSums(ann_tab))
  ann_tab <- as.data.frame(ann_tab)
  ann_tab_pct = lapply(ann_tab[,], function(x) {
    x/ann_tab$Total})
  ann_tab_pct <- as.data.frame(ann_tab_pct)
  rownames(ann_tab_pct) <- rownames(ann_tab)
  write.csv(ann_tab_pct, file = paste0(output_file_path,sample_name,"_proportions_",ident_of_interest1,"_",ident_of_interest2,".csv"))
}

### change the table to plot in barplot 
create_table_cell_type_prop_for_plot_one_sample <- function(
    df_cell_type_prop_per_sample, 
    sample_name,
    n_cell_types_dim
){
  df <- df_cell_type_prop_per_sample[df_cell_type_prop_per_sample$X %in% sample_name,]
  df <- t(df)
  df <- as.data.frame(df)
  colnames(df) <- NULL
  df <- df[n_cell_types_dim,]
  df$cell_types <- rownames(df)
  colnames(df) <- c("proportion","cell_types")
  rownames(df) <- NULL
  df$sample <- sample_name
  return(df)
}

create_table_cell_type_prop_table_for_plot <- function(
    df_cell_type_prop_per_sample, 
    n_cell_types_dim
) {
  df_cell_type_prop_per_sample$Total <- NULL
  
  sample_names <- df_cell_type_prop_per_sample$X
  df <- NULL
  
  for (i in sample_names){
    df1 <- create_table_cell_type_prop_for_plot_one_sample(df_cell_type_prop_per_sample,i,n_cell_types_dim) 
    df <- rbind(df, df1)
  }
  
  df$proportion <- as.numeric(df$proportion)
  return(df)
}

### statistical analysis 
cell_type_prop_stats <- function(
    seurat_object,
    anno_column,
    control,
    condition,
    condition_column,
    log2FD_cutoff,
    output_path
    
){
  #https://github.com/rpolicastro/scProportionTest
  prop_test <- sc_utils(seurat_object)
  #permutation testing and bootstrapping 
  #first sample is the control, the second sample is the condition 
  #negative obs_logFD means that the cell type is depleted in the condition 
  #positive obs_logFD means that the cell type is enriched in the condition 
  #grey is non-significant, red is significant 
  prop_test <- permutation_test(prop_test, cluster_identity = anno_column,
                                sample_1 = control, sample_2 = condition,sample_identity = condition_column)
  p <- permutation_plot(prop_test,log2FD_threshold = log2(log2FD_cutoff))
  print(p)
  ggsave(file = output_path, width = 8, height = 2, plot = p)
}

### barplot per cell type between two conditions
barplot_cell_type_oi_2cond <- function(
    cell_type_prop_df,
    col_names_oi,
    cols_oi,
    output_path
){
  df <- cell_type_prop_df[,colnames(cell_type_prop_df) %in% c(col_names_oi)]
  colnames(df) <- c("cond","prop")
  p <- ggplot(df, aes(fill=cond, y=prop, x=cond)) + theme_classic() +
    geom_bar(position="stack", stat="identity" ) + 
    theme(axis.text = element_text(size = 10)) + 
    theme(axis.title= element_text(size = 10))  + 
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))  + 
    xlab("Sample") + ylab("Cell type proportion") + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values=  cols_oi)
  print(p)
  ggsave(output_path, width = 4, height = 5, plot = p)
}

##### Functions for DEG analysis 
### Function to do DGE analysis between two conditions and writing the output to a csv file 
DEG_to_csv_two_cond <- function(
    seurat_object,
    cond1,
    cond2,
    only.pos_condition,
    logfc_threshold,
    csv_file_directory
){
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize",
                                 scale.factor = 10000,
                                 margin = 1, assay = "RNA")
  DefaultAssay(seurat_object) <- "RNA"
  markers <- FindMarkers(object = seurat_object, ident.1 = cond1, ident.2 = cond2, only.pos = only.pos_condition, min.pct = 0.25, 
                         logfc.threshold = logfc_threshold,slot = "data")
  write.csv(markers, file = csv_file_directory)
}

### Function to plot DEG results in volcano and label genes of interest 
volcano_DGE_showing_goi <- function(
    DEG_list_path, 
    avg_log2FC_cutoff,
    p_val_adj_cutoff,
    condition1,
    condition2,
    genes_of_interest,
    colors_dots,
    x_line_intercept,
    output_path
){
  DEG_list <- read.csv(DEG_list_path)
  DEG_list <- DEG_list %>% mutate(
    Expression = case_when(avg_log2FC >= avg_log2FC_cutoff & p_val_adj <= p_val_adj_cutoff ~ condition1,
                           avg_log2FC <= -avg_log2FC_cutoff & p_val_adj <= p_val_adj_cutoff ~ condition2,
                           TRUE ~ "Non sig."))
  
  rownames(DEG_list) <- DEG_list$X
  
  DEG_list$genelabels <- ifelse(rownames(DEG_list) %in% genes_of_interest,TRUE,FALSE)
  
  p <- ggplot(DEG_list, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
    geom_point(size=8, aes(color = Expression)) +
    geom_text_repel(aes(label=ifelse(DEG_list$genelabels, rownames(DEG_list),"")),size=5, max.overlaps = 1000) +
    xlab("logFC") + 
    ylab("-log10(p_val_adj)") + 
    scale_color_manual(values =c(colors_dots) ,guide = "none") + 
    theme(axis.title= element_text(size = 25)) + theme(axis.text = element_text(size = 30))  + 
    theme(plot.title = element_text(size = 25, face = "bold"))  + 
    theme(legend.title = element_text(size = 25), legend.text = element_text(size = 25)) + 
    geom_vline(xintercept = x_line_intercept, linetype = "dashed", color = "gray50") + #log Fold Change cutoff 
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray50") + #p value cutoff 
    theme_classic(base_size = 25) 
  print(p)
  ggsave(output_path, width = 12, height = 10, plot = p)
}

##### Functions for GO term analysis 
### Function to generate ranks per sample 
ranks_for_fGSEA <- function(
    DEG_file
){
  degs <- read.csv(file = DEG_file)
  degs <- na.omit(degs)
  degs$p_val_adj[degs$p_val_adj == 0] <- .Machine$double.xmin #replace 0 with lowest number if p-adjusted is 0
  print(sort(degs$p_val_adj))
  ranks <- degs %>% na.omit()%>% mutate(ranking=-log10(p_val_adj)/sign(avg_log2FC))
  print(sort(ranks$ranking)) 
  ranks <- ranks$ranking
  names(ranks) <- degs$X
  head(ranks, 10)
  #remove infinit numbers 
  ranks <- ranks[!is.infinite(ranks)]
}

### Function to run GSEA of KEGG, BP, Reactome and Hallmarks using fGSEA and saving the combined results in a csv file 
fGSEA_to_csv <- function(
    ranks_vector,
    output_file
){
  GO_BP <- fgsea(pathways = BP, 
                 stats = ranks_vector,
                 minSize=10,
                 maxSize=500,
                 nperm=1000000)
  BP_filtered <- GO_BP %>% filter(abs(NES)>0 & padj<0.05)  
  
  GO_Hallmarks <- fgsea(pathways = Hallmarks, 
                        stats = ranks_vector,
                        minSize=10,
                        maxSize=500,
                        nperm=1000000)
  Hallmarks_filtered <- GO_Hallmarks %>% filter(abs(NES)>0 & padj<0.05)  
  
  GO_KEGG <- fgsea(pathways = KEGG, 
                   stats = ranks_vector,
                   minSize=10,
                   maxSize=500,
                   nperm=1000000)
  KEGG_filtered <- GO_KEGG %>% filter(abs(NES)>0 & padj<0.05)  
  
  ## wirte output 
  BP_filtered_df <- as.data.frame(BP_filtered)
  Hallmarks_filtered_df <- as.data.frame(Hallmarks_filtered)
  KEGG_filtered_df <- as.data.frame(KEGG_filtered)
  
  all_filtered <- rbind(BP_filtered_df, Hallmarks_filtered_df)
  all_filtered <- rbind(all_filtered, KEGG_filtered_df)
  
  all_filtered <- apply(all_filtered,2,as.character)
  write.csv(all_filtered, file = output_file)
}

##### Functions to analyse MHC-I and antigen processing score between two conditions 
### mouse data 
MHCI_antigen_processing_vln_2_cond_mm <- function(
    seurat_object,
    column_oi1,
    ident_to_subset,
    condition_column,
    first_identity_for_stats, 
    second_identity_for_stats,
    colors_to_plot,
    y_limits,
    output_path_file
){
  Idents(seurat_object) <- column_oi1
  sub <- subset(seurat_object, idents = ident_to_subset)
  
  Antigen_list <- list(c("H2-D1",  "H2-Q7",  "H2-K1",  "H2-T23", "H2-Q4",  "Tap1",   "Tapbp",  "B2m"  , 
                         "Psmb8" , "Psme1",  "Psmb9",  "Calr" ,  "Psmb10", "Ncf1" ,  "Fcer1g"))
  
  sub <-AddModuleScore(sub, features= Antigen_list,name = "AntigenMHCI")
  names(x = sub[[]])
  
  #test
  Idents(sub) <- condition_column
  first <- subset(sub, idents = first_identity_for_stats)
  second <- subset(sub, idents = second_identity_for_stats)
  p_val <- wilcox.test(first$AntigenMHCI1, second$AntigenMHCI1, alternative = "two.sided")
  p_val <- p_val$p.value
  
  p <- VlnPlot(sub, features= "AntigenMHCI1", group.by = condition_column, pt.size = 0, cols = colors_to_plot) +  
    theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
    labs(title = paste0("cluster", ident_to_subset), y = " MHCI angigen pres and proc score", x="") + theme(legend.position="right") +  
    annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
    scale_y_continuous(limits = y_limits)
  print(p)
  ggsave(output_path_file, width = 8, height = 8, plot = p)
}

### human data 
MHCI_antigen_processing_vln_2_cond_hs <- function(
    seurat_object,
    column_oi1,
    ident_to_subset,
    condition_column,
    first_identity_for_stats, 
    second_identity_for_stats,
    colors_to_plot,
    y_limits,
    output_path_file
){
  Idents(seurat_object) <- column_oi1
  sub <- subset(seurat_object, idents = ident_to_subset)
  
  Antigen_list <- list(c("TAP1", "TAPBP", "B2M", "PSMB8", "PSMB1", "PSMB9", "CALR", "PSMB10", "NCF1", "FCER1G", "HLA-A", "HLA-C", "HLA-B"))
  
  sub <-AddModuleScore(sub, features= Antigen_list,name = "AntigenMHCI")
  names(x = sub[[]])
  
  #test
  Idents(sub) <- condition_column
  first <- subset(sub, idents = first_identity_for_stats)
  second <- subset(sub, idents = second_identity_for_stats)
  p_val <- wilcox.test(first$AntigenMHCI1, second$AntigenMHCI1, alternative = "two.sided")
  p_val <- p_val$p.value
  
  p <- VlnPlot(sub, features= "AntigenMHCI1", group.by = condition_column, pt.size = 0, cols = colors_to_plot) +  
    theme_classic() + theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
    labs(title = paste0("cluster", ident_to_subset), y = " MHCI angigen pres and proc score", x="") + theme(legend.position="right") +  
    annotate("text",x=1,y=0.5, label = paste0("p_val wilcox = ",p_val)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black") + 
    scale_y_continuous(limits = y_limits)
  print(p)
  ggsave(output_path_file, width = 8, height = 8, plot = p)
}

##### Functions for CellPhoneDB 
### Function for generation of input files from mouse data 
#first gene symbols have to be converted from mouse to human because CellPhoneDB database only contains human L-R interactions   
#Conversion of mouse to human genes adapted partially from: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md
#and: https://www.cellphonedb.org/faq-and-troubleshooting
#outputs are two text tiles: gene counts and cell annotations (meta)
Input_files_CellPhoneDB_generation_mm <- function(
    seurat_object,
    annotation_column,
    sample_name,
    ouput_file_path
){
  ### load human and mouse ensemble symbols
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  ###generating counts file 
  #take raw data and normalize it
  count_raw_meta <- GetAssayData(object = seurat_object, slot = "counts")[,colnames(x = seurat_object)]
  count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(count_norm_meta) , mart = mouse, 
                   attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
  print(head(genesV2))
  matrixA <- count_norm_meta[match(genesV2$MGI.symbol,rownames(count_norm_meta),nomatch=T),]
  matrixB <- matrixA
  matrixB$gene <- genesV2$Gene.stable.ID
  rownames(matrixA) <- matrixB$gene
  #save count matrix as text file 
  write.table(matrixA, paste0(ouput_file_path,sample_name,"_count.txt"), sep='\t', quote=F)
  ###generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(seurat_object@meta.data), seurat_object@meta.data[,annotation_column, drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, paste0(ouput_file_path,sample_name,"_meta.txt"), sep='\t', quote=F, row.names=F)
}

### Function for generation of CellPhoneDB input files from human data 
Input_files_CellPhoneDB_generation_hs <- function(
    seurat_object,
    annotation_column,
    sample_name,
    ouput_file_path
){
  ###generating counts file 
  #take raw data and normalize it
  count_raw_meta <- GetAssayData(object = seurat_object, layer = "counts")[,colnames(x = seurat_object)]
  count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
  write.table(count_norm_meta, paste0(ouput_file_path,sample_name,"_count.txt"), sep='\t', quote=F)
  ###generating meta file based on cell type annotation of Seurat object 
  meta_data_meta <- cbind(rownames(seurat_object@meta.data), seurat_object@meta.data[,annotation_column, drop=F])  
  #save meta file as text file 
  write.table(meta_data_meta, paste0(ouput_file_path,sample_name,"_meta.txt"), sep='\t', quote=F, row.names=F)
}

