library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(data.table)

sharing = readr::read_tsv("mash_sharing.tsv")

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

# friendly label
ontology_map = ontology_map %>% dplyr::mutate(label = paste(study, ontology_tissue, sep=" "))

# lead datasets
leads=c("TwinsUK.blood", "BLUEPRINT_PE.T-cell", "GEUVADIS.LCL", 
        "BLUEPRINT_SE.neutrophil", "BLUEPRINT_SE.monocyte", "Alasoo_2018.macrophage_naive", 
        "ROSMAP.brain_naive", "HipSci.iPSC")

plot_distribution <- function(similarity_matrix, y_axis="Mash sharing", measure="sharing"){
  df_sharing = pivot_longer(similarity_matrix, cols=c(-dataset), names_to="dataset2")
  # remove diagonal values
  df_sharing = df_sharing %>% filter(value < 1)
  # filter similarities for lead datasets
  df_sharing = filter(df_sharing, dataset %in% leads)
  df_sharing = df_sharing %>% left_join(ontology_map[c("tissue_fct", "study_qtlgroup")], by=c("dataset2"="study_qtlgroup"))
  # add friendly label to lead datasets
  df_sharing = df_sharing %>% left_join(ontology_map[c("label", "study_qtlgroup")], by=c("dataset"="study_qtlgroup"))
  
  df_sharing = dplyr::rename(df_sharing, !!measure:=value)
  
  # or ggplot(df_sharing, aes(dataset, get(measure), colour=tissue_fct)) for friendly label of lead group
  plt <- ggplot(df_sharing, aes(dataset, get(measure), colour=tissue_fct)) +
    geom_jitter(width = 0.2) +
    xlab("Dataset") +
    ylab(y_axis) +
    scale_colour_manual(name = "group",
                        values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#fed976","#f781bf","#999999"))  +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(plt)
}


plt = plot_distribution(sharing)

ggsave("sharing_distribution.pdf", plt, width = 8, height = 8)

