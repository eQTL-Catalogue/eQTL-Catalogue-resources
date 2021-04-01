library(tidyverse)
library(ggplot2)

sample_size = readr::read_tsv("dataset_sample_sizes.tsv")
sample_size = sample_size %>% mutate(study = ifelse(study %in% c("BLUEPRINT_PE","BLUEPRINT_SE"), "BLUEPRINT", study))

ontology = readr::read_tsv("../../../ontology_mappings/tissue_ontology_mapping.tsv")
friendly_names = readr::read_tsv("../../../ontology_mappings/friendly_names.tsv")

ontology = ontology %>% dplyr::left_join(friendly_names[c("tissue_ontology_id", "tissue_label")], by="tissue_ontology_id")
ontology = ontology %>% dplyr::arrange(tissue_label, study)

studies = dplyr::left_join(sample_size, ontology)

rnaseq_studies = dplyr::filter(studies, type == "RNA-seq")
microarray_studies = dplyr::filter(studies, type == "microarray", study != "Kolberg_2020")

# set GTEx to the bottom of plot stacks
levels = c(sort(setdiff(unique(rnaseq_studies$study), "GTEx")), "GTEx")
rnaseq_studies = rnaseq_studies %>% 
  mutate(study = factor(study, levels))

# assign default colours to studies and dark grey to GTEx 
gg_color_hue <- function(n) {
  # https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

adj_names = sort(setdiff(unique(rnaseq_studies$study), "GTEx"))
values = gg_color_hue(length(adj_names))
names(values) = adj_names
clr_values = c(values, c(GTEx="darkgrey"))


draw_plot = function(studies, legend_rows = 4, legend_x_pos = 0.5, legen_y_position=0.8){
  sample_sizes = studies %>% 
    dplyr::group_by(tissue_label) %>% 
    dplyr::mutate(total_sample_size=sum(dataset_sample_size))
  
  plt = ggplot(sample_sizes, aes(x = reorder(tissue_label, -total_sample_size), dataset_sample_size, fill=study)) +
    geom_col() + 
    guides(fill=guide_legend(nrow=legend_rows,byrow=TRUE))+
    xlab("Cell type or tissue") + 
    ylab("Sample size") +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
          panel.grid = element_blank(),
          legend.position=c(legend_x_pos, legen_y_position),
          legend.background = element_rect(colour="lightgrey", 
                                           size=0.5, linetype="solid"))
  
  if("GTEx" %in% sample_sizes$study){
    plt = plt + scale_fill_manual(values=clr_values)
  }
  return(plt)
}

rnaseq_plt = draw_plot(rnaseq_studies, legen_y_position = 0.75)
ggsave("rnaseq_sample_size.pdf", rnaseq_plt, width = 10, height = 5)

microarr_plt = draw_plot(microarray_studies, legend_rows = 3, legend_x_pos = 0.76, legen_y_position=0.76)
ggsave("microarr_sample_size.pdf", microarr_plt, width = 6, height = 3.7)
