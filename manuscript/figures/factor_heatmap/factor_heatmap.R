library(tidyverse)
library(data.table)

# load factor matrix
FactorM = read.table("K30_a1900_l11100_Run4.txt")

# sorting and naming of studies
###############################
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
# place immune t-cells together
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(qtl_group %like% "anti", "ANTI", ontology_tissue), 
                                              ontology_tissue = ifelse(qtl_group %like% "anti", paste(ontology_tissue, "(anti-CD3-CD28)"), ontology_tissue))
# place quach monocyte and nedelec macrophage together because they load on one factor
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(study_qtlgroup == "Quach_2016.monocyte_naive", "z-monocyte", dummy))
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(study_qtlgroup == "Nedelec_2016.macrophage_naive", "a-macrophage", dummy))
# sort studies
ontology_map = ontology_map %>% dplyr::arrange(tissue_fct, dummy, study)
ontology_map = ontology_map %>% dplyr::mutate(heatmap_label = paste(study, ontology_tissue, sep=" "))
ontology_map = ontology_map %>% dplyr::filter(study_qtlgroup %in% rownames(FactorM))

labels = data.frame(ontology = ontology_map$heatmap_label)
rownames(labels) = ontology_map$study_qtlgroup

# sort tissues
FactorM = FactorM[ontology_map$study_qtlgroup, ]
# remove stimulated datasets
datasets = rownames(FactorM)
datasets = datasets[!(datasets %like% "IFNg" | 
                        datasets %like% "Salmonella" | 
                        datasets %like% "Listeria" | 
                        datasets %like% "IAV" | 
                        datasets %like% "Pam3CSK4" | 
                        datasets %like% "LPS" | 
                        datasets %like% "R848")]
FactorM = FactorM[datasets,]

# name factors
factor_names = tibble(factor = paste0("Factor", c(1:16)),
                      name = c("Universal", "iPSC", "Skin", "Blood",
                               "LCL", "T-cell", "Monocyte & Macrophage", "Fat",
                               "T-cell immune response", "Brain", "Macrophage", "Monocyte naive",
                               "BLUEPRINT T-cell", "ROSMAP Brain", "Muscle", "Neutrophil" ))
colnames(FactorM) = factor_names$name

colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#fed976","#f781bf","#999999")
names(colors) = c("blood","lymphocyte","LCL","neutrophil","monocyte","macrophage","brain","iPSC","other")

# column annotation for heatmap
col_annot = data.frame(group = ontology_map$tissue_fct)
rownames(col_annot) = ontology_map$study_qtlgroup

pheatmap::pheatmap(t(FactorM), fontsize=8, cluster_rows = F, cluster_cols = F, cellwidth=10, cellheight = 10,
                   angle_col = 45, labels_col = labels[rownames(FactorM),], 
                   annotation_legend = T, annotation_col = col_annot, annotation_colors=list(group=colors),
                   height = 6, width = 17, filename="factor_heatmap.pdf")


# remove datasets with small loadings across non-universal factors
matrix = FactorM[rowSums(FactorM[,2:16]) > 3,]

pheatmap::pheatmap(t(matrix), fontsize=8, cluster_rows = F, cluster_cols = F, cellwidth=10, cellheight = 10,
                   angle_col = 45, labels_col = labels[rownames(matrix),], 
                   annotation_legend = T, annotation_col = col_annot, annotation_colors=list(group=colors), 
                   height = 6, width = 13, filename="factor_heatmap_compact.pdf")

