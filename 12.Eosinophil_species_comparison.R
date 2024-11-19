########## This code compares eosinophils between human and mouse ##########

##### link to libraries and functions
source("~/Projects/Eosinophils_in_late_stage_CRC/1.Packages_and_functions.R")

##### load R object 
human_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/human_eosinophils_integrated_annotated.rds")
mouse_eos <- readRDS("/data/khandl/Eosinophils_in_CRC/seurat_objects/mouse_eosinophils_integrated_annotated.rds")

##### compare tissues in human data with conditions in mouse data 
### generate dataframes with average expression per condition/tissue 
average_expression_human <- AverageExpression(human_eos, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "tissue")
average_expression_human_df <- as.data.frame(average_expression_human)

average_expression_mouse <- AverageExpression(mouse_eos, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "condition")
average_expression_mouse_df <- as.data.frame(average_expression_mouse)

### Venn diagrams 
## blood 
# remove 0 counts 
average_expression_human <- average_expression_human_df[average_expression_human_df$RNA.blood.healhty >0 ,]
average_expression_mouse <- average_expression_mouse_df[average_expression_mouse_df$RNA.blood.wt >0 ,]

# remove where no orthologues are found 
#human data 
mouse_orthologs_from_human_data <- convert_orthologs(gene_df = average_expression_human,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "human",
                                                     output_species = "mouse",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_human_remove_no_orthologues <- average_expression_human[!rownames(average_expression_human) %in% rownames(mouse_orthologs_from_human_data),]

#mouse data 
human_orthologs_from_mouse_data <- convert_orthologs(gene_df = average_expression_mouse,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "mouse",
                                                     output_species = "human",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_mouse_remove_no_orthologues <- average_expression_mouse[!rownames(average_expression_mouse) %in% rownames(human_orthologs_from_mouse_data),]

x <- list("M" = rownames(human_orthologs_from_mouse_data), "H" =rownames(average_expression_human_remove_no_orthologues))
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

## tumor 
# remove 0 counts 
average_expression_human <- average_expression_human_df[average_expression_human_df$RNA.tumor >0 ,]
average_expression_mouse <- average_expression_mouse_df[average_expression_mouse_df$RNA.tumor.wt >0 ,]

# remove where no orthologues are found 
#human data 
mouse_orthologs_from_human_data <- convert_orthologs(gene_df = average_expression_human,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "human",
                                                     output_species = "mouse",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_human_remove_no_orthologues <- average_expression_human[!rownames(average_expression_human) %in% rownames(mouse_orthologs_from_human_data),]

#mouse data 
human_orthologs_from_mouse_data <- convert_orthologs(gene_df = average_expression_mouse,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "mouse",
                                                     output_species = "human",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_mouse_remove_no_orthologues <- average_expression_mouse[!rownames(average_expression_mouse) %in% rownames(human_orthologs_from_mouse_data),]

x <- list("M" = rownames(human_orthologs_from_mouse_data), "H" =rownames(average_expression_human_remove_no_orthologues))
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

## control 
# remove 0 counts 
average_expression_human <- average_expression_human_df[average_expression_human_df$RNA.tissue.ctrl >0 ,]
average_expression_mouse <- average_expression_mouse_df[average_expression_mouse_df$RNA.adult.colon.wt >0 ,]

# remove where no orthologues are found 
#human data 
mouse_orthologs_from_human_data <- convert_orthologs(gene_df = average_expression_human,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "human",
                                                     output_species = "mouse",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_human_remove_no_orthologues <- average_expression_human[!rownames(average_expression_human) %in% rownames(mouse_orthologs_from_human_data),]

#mouse data 
human_orthologs_from_mouse_data <- convert_orthologs(gene_df = average_expression_mouse,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "mouse",
                                                     output_species = "human",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_mouse_remove_no_orthologues <- average_expression_mouse[!rownames(average_expression_mouse) %in% rownames(human_orthologs_from_mouse_data),]

x <- list("M" = rownames(human_orthologs_from_mouse_data), "H" =rownames(average_expression_human_remove_no_orthologues))
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

##### compare A eos and B eos from different species 
### generate dataframes with average expression per condition/tissue 
average_expression_human <- AverageExpression(human_eos, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "annotation")
average_expression_human_df <- as.data.frame(average_expression_human)

average_expression_mouse <- AverageExpression(mouse_eos, return.seurat = FALSE, normalization.method = "LogNormalize",assays = "RNA", group.by = "annotation")
average_expression_mouse_df <- as.data.frame(average_expression_mouse)

### Venn diagrams 
## A-eos  
# remove 0 counts 
average_expression_human <- average_expression_human_df[average_expression_human_df$RNA.A.eos >0 ,]
average_expression_mouse <- average_expression_mouse_df[average_expression_mouse_df$RNA.A.eos >0 ,]

# remove where no orthologues are found 
#human data 
mouse_orthologs_from_human_data <- convert_orthologs(gene_df = average_expression_human,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "human",
                                                     output_species = "mouse",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_human_remove_no_orthologues <- average_expression_human[!rownames(average_expression_human) %in% rownames(mouse_orthologs_from_human_data),]

#mouse data 
human_orthologs_from_mouse_data <- convert_orthologs(gene_df = average_expression_mouse,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "mouse",
                                                     output_species = "human",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_mouse_remove_no_orthologues <- average_expression_mouse[!rownames(average_expression_mouse) %in% rownames(human_orthologs_from_mouse_data),]

x <- list("M" = rownames(human_orthologs_from_mouse_data), "H" =rownames(average_expression_human_remove_no_orthologues))
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

## B-eos  
# remove 0 counts 
average_expression_human <- average_expression_human_df[average_expression_human_df$RNA.B.eos >0 ,]
average_expression_mouse <- average_expression_mouse_df[average_expression_mouse_df$RNA.B.eos >0 ,]

# remove where no orthologues are found 
#human data 
mouse_orthologs_from_human_data <- convert_orthologs(gene_df = average_expression_human,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "human",
                                                     output_species = "mouse",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_human_remove_no_orthologues <- average_expression_human[!rownames(average_expression_human) %in% rownames(mouse_orthologs_from_human_data),]

#mouse data 
human_orthologs_from_mouse_data <- convert_orthologs(gene_df = average_expression_mouse,
                                                     gene_input = "rownames", 
                                                     gene_output = "rownames", 
                                                     input_species = "mouse",
                                                     output_species = "human",
                                                     non121_strategy = "drop_both_species",
                                                     method = "gprofiler") 
#remove genes where no orhtologs are found 
average_expression_mouse_remove_no_orthologues <- average_expression_mouse[!rownames(average_expression_mouse) %in% rownames(human_orthologs_from_mouse_data),]

x <- list("M" = rownames(human_orthologs_from_mouse_data), "H" =rownames(average_expression_human_remove_no_orthologues))
ggVennDiagram(x) + theme(plot.title = element_text(size = 25, face = "bold")) 

