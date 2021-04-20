library("tidyverse")
library("ggplot2")
library("data.table")

# path_files <- list.files("../../temp_files/eQTL_sharing/new_with_GTExV8/txrevise_no_graph_mash/mash/", full.names = T)
# path_files <- "../../temp_files/eQTL_sharing/new_with_GTExV8/exon_no_graph_mash//"
path_files <- "../../temp_files/eQTL_sharing/new_with_GTExV8/ge_graph_30_thresh///"
# path_files <- "../../temp_files/eQTL_sharing/new_with_GTExV8/ma_graph_30_thresh//"
# path_files <- "/Users/kerimov/Work/temp_files/eQTL_sharing/new_with_GTExV8/txrev_no_graph_all_together/"
# path_files <- "/Users/kerimov/Work/temp_files/eQTL_sharing/new_with_GTExV8/ge_graph_30_minpval/"
for (path_file in path_files) {
  load(paste0(path_file,"/sharing.R"))
  # sharing = read_tsv("../../temp_files/eQTL_sharing/results/results_min_pvalue/sharing.tsv")
  # sharing_matrix = as.matrix(sharing %>% select(-dataset))
  sharing_matrix = as.matrix(sharing)
  rownames(sharing_matrix) = gsub(pattern = "GTExV8", replacement = "GTEx", x = colnames(sharing_matrix))
  colnames(sharing_matrix) <- gsub(pattern = "GTExV8", replacement = "GTEx", x = colnames(sharing_matrix))
  fit <- MASS::isoMDS(1-sharing_matrix, k=2)
  
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
  ontology_map = ontology_map %>% dplyr::filter(study_qtlgroup %in% rownames(sharing_matrix))
  
  plot_coords = function(coords){
    coords = tibble(x = coords$points[,1], y=coords$points[,2], study_qtlgroup=rownames(coords$points))
    coords = left_join(coords, ontology_map[c("study_qtlgroup", "tissue_fct")])
    coords$study = ifelse(grepl("GTEx", coords$study_qtlgroup), "GTEx", NA)
    coords$study = ifelse(grepl("BLUEPRINT", coords$study_qtlgroup), "BLUEPRINT", coords$study)
    coords$study = factor(coords$study, levels=c("GTEx", "BLUEPRINT", NA))
    plt = ggplot2::ggplot(coords, aes(x, y, grp=study_qtlgroup)) +
      ggplot2::geom_point(aes(colour=tissue_fct), size=3, show.legend = F) +
      geom_point(aes(shape=study, size=study), data=coords[!is.na(coords$study),], fill=NA) +
      scale_shape_manual(values=c(21,24))+
      scale_size_manual(values=c(3,5))+
      scale_colour_manual(name="group", values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#fed976","#f781bf","#999999"))  +  
      theme_light() + 
      theme(panel.grid = element_blank()) + 
      ggplot2::labs(x="MDS Coordinate 1", y="MDS Coordinate 2")
    return(plt)
  }
  
  fit$points <- -fit$points
  
  plt = plot_coords(fit)
  ggsave(paste0("mash_mds_inverted.pdf"), plt, width = 5, height = 3.1)
  
  pltly <- plotly::ggplotly(plt)
  htmlwidgets::saveWidget(widget = plotly::as_widget(pltly),
                          file = file.path("mash_mds_inverted.html"),
                          libdir = "dependencies")
}
