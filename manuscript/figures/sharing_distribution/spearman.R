library(tidyverse)
library(data.table)
library(RColorBrewer)

path_files <- list.files("../../temp_files/eQTL_sharing/new_with_GTExV8/txrevise_no_graph_mash/similarity//", full.names = T)
path_files <- c(path_files,"../../temp_files/eQTL_sharing/new_with_GTExV8/tx_no_graph_mash/")
path_files <- c(path_files,"../../temp_files/eQTL_sharing/new_with_GTExV8/ge_graph_15_thresh/")
path_files <- c(path_files,"../../temp_files/eQTL_sharing/new_with_GTExV8/exon_no_graph_mash/")
path_files <- c(path_files,"../../temp_files/eQTL_sharing/new_with_GTExV8/ge_graph_30_thresh/")
path_files <- c(path_files,"/Users/kerimov/Work/temp_files/eQTL_sharing/new_with_GTExV8/ge_no_graph_mash/")
path_files <- c("/Users/kerimov/Work/temp_files/eQTL_sharing/new_with_GTExV8/txrev_no_graph_all_together//")
path_files <- c("/Users/kerimov/Work/temp_files/eQTL_sharing/new_with_GTExV8/ge_graph_30_minpval/")

for (path_file in path_files) {
  correlation_matrix = read_delim(paste0(path_file, "/spearman_cor_na_to_zero.txt"), delim = " ", skip = 1, col_names = F)
  correlation_matrix$X1 <- gsub(pattern = "GTExV8", replacement = "GTEx", x = correlation_matrix$X1)
  colnames(correlation_matrix) <- c("dataset",correlation_matrix$X1)
  # correlation_matrix = read_tsv("../../temp_files/eQTL_sharing/results/results_random/spearman_cor_na_to_zero.txt")
  correlation_matrix = as.matrix(correlation_matrix %>% select(-dataset))
  rownames(correlation_matrix) = colnames(correlation_matrix)
  
  # naming of studies
  ###############################
  ontology_map = readr::read_tsv("../../../ontology_mappings/tissue_ontology_mapping.tsv")
  friendly_names = readr::read_tsv("../../../ontology_mappings/friendly_names.tsv") %>%
    dplyr::select(tissue_ontology_id, tissue_label)
  ontology_map$study[c(1,2,3)] = c("BLUEPRINT_SE", "BLUEPRINT_SE", "BLUEPRINT_PE")
  ontology_map <- ontology_map %>% dplyr::mutate(study_qtlgroup = paste0(study, ".", qtl_group)) %>%
    dplyr::left_join(friendly_names) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "brain", "brain", "other")) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "LCL", "LCL", sample_class)) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "monocyte", "monocyte", sample_class)) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "macrophage", "macrophage", sample_class)) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "blood", "blood", sample_class)) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "neutrophil", "neutrophil", sample_class)) %>%
    dplyr::mutate(sample_class = ifelse(tissue_ontology_id %in% c("CL_0000236","CL_0002677","CL_0002678","CL_0000624","CL_0000625","CL_0000623","CL_0000899","CL_0000546","CL_0000545","CL_0000899","CL_0002038","CL_0000084"), "lymphocyte", sample_class)) %>%
    dplyr::mutate(sample_class = ifelse(tissue_label %like% "iPSC", "iPSC", sample_class))
  
  fct_levels = c("blood","lymphocyte","LCL","neutrophil","monocyte","macrophage","brain","iPSC","other")
  ontology_map = dplyr::mutate(ontology_map, tissue_fct = factor(sample_class, levels = fct_levels))
  ontology_map = ontology_map %>% dplyr::mutate(heatmap_label = paste(study, tissue_label, sep=" "))
  ontology_map = ontology_map %>% dplyr::filter(study_qtlgroup %in% rownames(correlation_matrix))
  
  labels = data.frame(ontology = ontology_map$heatmap_label)
  rownames(labels) = ontology_map$study_qtlgroup
  
  colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#fed976","#f781bf","#999999")
  names(colors) = c("blood","lymphocyte","LCL","neutrophil","monocyte","macrophage","brain","iPSC","other")
  
  # column annotation for heatmap
  col_annot = data.frame(group = ontology_map$tissue_fct)
  rownames(col_annot) = ontology_map$study_qtlgroup
  
  # rownames(correlation_matrix) = labels[rownames(correlation_matrix),]
  
  pheatmap::pheatmap(correlation_matrix, fontsize=8, cluster_rows = T, cluster_cols = T, cellwidth = 10, cellheight = 8,
                     color=colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100), border_color=NA,
                     angle_col = 45,   
                     annotation_legend = T, annotation_col = col_annot, annotation_colors=list(group=colors), 
                     filename = paste0(path_file,"/spearman_heatmap.pdf"), width = 20, height = 14)
  
  # pheatmap::pheatmap(correlation_matrix, fontsize=8, cluster_rows = T, cluster_cols = T, cellwidth = 10, cellheight = 8,
  #                    color=colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100), border_color=NA,
  #                    angle_col = 45, cutree_rows = 10,
  #                    annotation_legend = T, annotation_col = col_annot, annotation_colors=list(group=colors), 
  #                    filename = paste0("mash_heatmap.pdf"), width = 20, height = 14)
}  
  
