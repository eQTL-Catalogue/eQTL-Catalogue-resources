library(tidyverse)
library(ggplot2)

sample_size = readr::read_tsv("dataset_sample_sizes.tsv")

ontology = readr::read_tsv("../../../ontology_mappings/tissue_ontology_mapping.tsv")
friendly_names = readr::read_tsv("../../../ontology_mappings/friendly_names.tsv")

ontology = ontology %>% dplyr::left_join(friendly_names[c("ontology_term", "ontology_tissue")], by="ontology_term")
ontology = ontology %>% dplyr::arrange(ontology_tissue, study)

studies = dplyr::left_join(sample_size, ontology)

rnaseq_studies = dplyr::filter(studies, type == "RNA-seq")
microarray_studies = dplyr::filter(studies, type == "microarray", study != "Kolberg_2020")

draw_plot = function(studies, legend_rows = 4, legen_y_position=0.8){
  sample_sizes = studies %>% 
    dplyr::group_by(ontology_tissue) %>% 
    dplyr::mutate(total_sample_size=sum(dataset_sample_size))
  
  ggplot(sample_sizes, aes(x = reorder(ontology_tissue, -total_sample_size), dataset_sample_size, fill=study)) +
    geom_col() + 
    guides(fill=guide_legend(nrow=legend_rows,byrow=TRUE))+
    xlab("Tissue (naive)") + 
    ylab("Sample size") +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
          panel.grid = element_blank(),
          legend.position=c(0.5, legen_y_position),
          legend.background = element_rect(colour="lightgrey", 
                                           size=0.5, linetype="solid"))
}
#6x12inch
rnaseq_plt = draw_plot(rnaseq_studies)
ggsave("rnaseq_sample_size.pdf", rnaseq_plt, width = 12, height = 6)
#6x9
microarr_plt = draw_plot(microarray_studies, legend_rows = 1, legen_y_position=0.9)
ggsave("microarr_sample_size.pdf", microarr_plt, width = 9, height = 6)
