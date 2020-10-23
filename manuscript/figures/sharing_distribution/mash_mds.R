library(tidyverse)
library(ggplot2)

sharing = read_tsv("mash_sharing.tsv")
sharing_matrix = as.matrix(sharing %>% select(-dataset))
rownames(sharing_matrix) = colnames(sharing_matrix)
fit <- MASS::isoMDS(1-sharing_matrix, k=2)

ontology_map = readr::read_tsv("../../../ontology_mappings/tissue_ontology_mapping.tsv")
friendly_names = readr::read_tsv("../../../ontology_mappings/friendly_names.tsv") %>%
  dplyr::select(ontology_term, ontology_tissue)
ontology_map$study[c(1,2,3)] = c("BLUEPRINT_SE", "BLUEPRINT_SE", "BLUEPRINT_PE")
ontology_map <- ontology_map %>% dplyr::mutate(study_qtlgroup = paste0(study, ".", qtl_group)) %>%
  dplyr::left_join(friendly_names) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "brain", "brain", "other")) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "LCL", "LCL", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "monocyte", "monocyte", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "macrophage", "macrophage", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "blood", "blood", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "neutrophil", "neutrophil", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_term %in% c("CL_0000236","CL_0002677","CL_0002678","CL_0000624","CL_0000625","CL_0000623","CL_0000899","CL_0000546","CL_0000545","CL_0000899","CL_0002038","CL_0000084"), "lymphocyte", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "iPSC", "iPSC", sample_class))

fct_levels = c("blood","lymphocyte","LCL","neutrophil","monocyte","macrophage","brain","iPSC","other")
ontology_map = dplyr::mutate(ontology_map, tissue_fct = factor(sample_class, levels = fct_levels))
ontology_map = ontology_map %>% dplyr::mutate(heatmap_label = paste(study, ontology_tissue, sep=" "))
ontology_map = ontology_map %>% dplyr::filter(study_qtlgroup %in% rownames(sharing_matrix))

plot_coords = function(coords){
  coords = tibble(x = coords$points[,1], y=coords$points[,2], study_qtlgroup=rownames(coords$points))
  coords = left_join(coords, ontology_map[c("study_qtlgroup", "tissue_fct")])
  coords
  plt = ggplot2::ggplot(coords, aes(x, y, colour=tissue_fct, grp=study_qtlgroup)) +
    ggplot2::geom_point() +
    scale_colour_manual(name = "group",
                        values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#fed976","#f781bf","#999999"))  +
    theme_light() + 
    theme(panel.grid = element_blank())
  return(plt)
}

plt = plot_coords(fit)
ggsave("mash_mds.pdf", plt, width = 8, height = 6)
