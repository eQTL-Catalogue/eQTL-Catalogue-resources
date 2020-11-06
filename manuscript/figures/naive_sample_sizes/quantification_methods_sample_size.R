library(tidyverse)
library(ggplot2)
# library(plotly)

#################################

# sumstat_path = "/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/"
# suffix = ".permuted.txt.gz" 
# files <- list.files(sumstat_path, pattern=suffix)
# 
# datasets = lapply(files, function(file){
#   data = readr::read_tsv(file.path(sumstat_path, file))
#   data$p_adj = p.adjust(data$p_beta, method = "fdr")
#   n = nrow(dplyr::filter(data, p_adj <= 0.05))
#   return(list(dataset=sub(".permuted.txt.gz", "", file), n_significant=n))
# })
# 
# datasets = dplyr::bind_rows(datasets)
# readr::write_tsv(datasets, "dataset_significant_qtls.tsv")

#################################
# sample sizes
# dir_path = "../SampleArcheology/studies/cleaned"
# files = list.files(dir_path, full.names = F)
# files = files[!files %in% c("GEUVADIS.tsv", "BLUEPRINT.tsv")]
# tbls = lapply(files, function(file){
#   tbl = read_tsv(file.path(dir_path, file), col_types = cols_only(rna_qc_passed="l", 
#                                                                   genotype_qc_passed="l", 
#                                                                   sample_id="c", 
#                                                                   qtl_group="c", 
#                                                                   cell_type="c", 
#                                                                   type="c", 
#                                                                   study="c"))
#   return(tbl)
# })
# studies = bind_rows(tbls)
# studies = studies %>% filter(rna_qc_passed == TRUE, genotype_qc_passed == TRUE)
# sample_size = studies %>% group_by(qtl_group, study) %>% summarise(dataset_sample_size=n())
# sample_size = sample_size %>% mutate(dataset=paste(study, qtl_group, sep="."))
# sample_size %>% write_tsv("sample_sizes.tsv")

#################################
sample_size = readr::read_tsv("sample_sizes.tsv")
datasets = readr::read_tsv("dataset_significant_qtls.tsv")
# split into datast and quantification method 
datasets = datasets %>% tidyr::separate(dataset, into=c("dataset", "term"), sep="_(?!.*_)")


datasets = datasets %>% left_join(sample_size)
datasets = datasets %>% mutate(type=if_else(term=="microarray", "Microarray", "RNA-seq")) %>% 
  mutate(term=if_else(term=="microarray", "ge", term))

datasets$term = factor(datasets$term, levels = c("ge", "exon", "tx", "txrev"),
                       labels = c("Gene expression", "Exon expression", "Transcript usage", "Event usage"))

datasets$type = factor(datasets$type, levels=c("RNA-seq", "Microarray"))
# 4x12  
plt = ggplot(datasets, aes(dataset_sample_size, n_significant, colour=type, dataset=dataset)) +
  geom_point() +
  facet_grid(. ~ term) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(), panel.grid = element_blank())+
  ylab("Number of QTLs") + 
  xlab("Sample size")

ggsave("quantification_methods_sample_size.pdf", plt, width=8, height = 2.8)

# ggplotly(plt)
