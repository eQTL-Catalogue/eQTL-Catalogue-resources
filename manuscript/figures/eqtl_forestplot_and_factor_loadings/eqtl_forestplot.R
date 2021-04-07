library(tidyverse)
library(ggplot2)
library(data.table)

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
# place immune t-cells together
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(qtl_group %like% "anti", "ANTI", tissue_label), 
                                              tissue_label = ifelse(qtl_group %like% "anti", paste(tissue_label, "(anti-CD3-CD28)"), tissue_label))
# place quach monocyte and nedelec macrophage together because they load on one factor
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(study_qtlgroup == "Quach_2016.monocyte_naive", "z-monocyte", dummy))
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(study_qtlgroup == "Nedelec_2016.macrophage_naive", "a-macrophage", dummy))
# sort studies
ontology_map = ontology_map %>% dplyr::arrange(tissue_fct, dummy, study)
ontology_map = ontology_map %>% dplyr::mutate(heatmap_label = paste(study, tissue_label, sep=" "))

# remove stimulated datasets
datasets = ontology_map$study_qtlgroup
datasets = !(datasets %like% "IFNg" | 
               datasets %like% "Salmonella" | 
               datasets %like% "Listeria" | 
               datasets %like% "IAV" | 
               datasets %like% "Pam3CSK4" | 
               datasets %like% "LPS" | 
               datasets %like% "R848")
ontology_map = ontology_map[datasets,]

# read effects data for the RBMS1 gene
effects = readr::read_tsv("rbms1_effects.tsv")
eqtl = effects %>% dplyr::filter(eqtl_id == "chr2_160430631_G_A.ENSG00000153250")

# function to extract beta, se and replace NA if needed
effects_to_matricies = function(effects, replace_na_with = FALSE){
  effects_matrix <- dplyr::select(effects, ends_with('.beta')) %>% 
    dplyr::rename_all(function(x) {sub(".beta", "", x)})
  errors_matrix <- dplyr::select(effects, ends_with('.se')) %>% 
    dplyr::rename_all(function(x){sub(".se", "", x)})
  
  effects_matrix <- as.matrix(effects_matrix)
  errors_matrix <- as.matrix(errors_matrix)
  
  missing_values <- which(is.na(effects_matrix), arr.ind=TRUE)
  if(replace_na_with == "mean"){
    effects_matrix[missing_values] <- rowMeans(effects_matrix, na.rm=TRUE)[missing_values[,1]]
    errors_matrix[missing_values] <- rowMeans(errors_matrix, na.rm=TRUE)[missing_values[,1]]
  }else if(replace_na_with == "zero"){
    effects_matrix[missing_values] <-  0
    errors_matrix[missing_values] <- 1
  }
  return(list(beta=effects_matrix, se=errors_matrix))
}

data = effects_to_matricies(eqtl)
effect = tibble(beta = data$beta[1,], qtl_group = names(data$beta[1,]), se=data$se[1,])
# calculate 95% confidence interval
conf.level = 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
effect = mutate(effect, interval = ci.value * se) %>% 
  mutate(qtl_group = gsub(pattern = "GTExV8", replacement = "GTEx", x = qtl_group))


colors = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#fed976","#f781bf","#999999")
names(colors) = c("blood","lymphocyte","LCL","neutrophil","monocyte","macrophage","brain","iPSC","other")

effect = dplyr::inner_join(ontology_map[c("heatmap_label", "study_qtlgroup", "tissue_fct")], effect, by = c("study_qtlgroup"="qtl_group"))
effect = effect %>% dplyr::mutate(study_qtlgroup=factor(study_qtlgroup, levels=study_qtlgroup))
plt = ggplot(effect, aes(x = study_qtlgroup, y = beta, ymin = beta - interval, ymax = beta + interval, colour=tissue_fct)) + 
  geom_point() + 
  geom_hline(yintercept=0, colour="grey") + 
  geom_errorbar(width = 0.1) + 
  scale_color_manual(values=colors, name="Group")+
  scale_x_discrete(labels = effect$heatmap_label) + 
  xlab("Dataset") + 
  ylab("Effect size") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(1, 1, 1, 3, "cm"),panel.grid = element_blank(),)+
  geom_hline(yintercept = 0)
ggsave("RBMS1_effects.pdf", plt, width = 10, height = 3.3)
