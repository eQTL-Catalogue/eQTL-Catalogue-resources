library("dplyr")

#Import mappings
friendly_names = readr::read_tsv("ontology_mappings/friendly_names.tsv") %>%
  dplyr::select(tissue_ontology_id, tissue_label)
mappings = readr::read_tsv("ontology_mappings/tissue_ontology_mapping.tsv") %>%
  dplyr::left_join(friendly_names, by = "tissue_ontology_id")

#Find those conditions that are not naive
condition_mapping = readr::read_tsv("ontology_mappings/cell_type_condition_mapping.tsv") %>%
  dplyr::select(study, qtl_group, condition_label) %>%
  dplyr::filter(condition_label != "naive")

#Add quants
quants = dplyr::tibble(quant_method = c("ge","exon","tx","txrev"))

#Specify microarray studies
microarray = dplyr::tibble(study = c("CEDAR", "Kasela_2017", 
                "Fairfax_2014", "Fairfax_2012", 
                "Naranbhai_2015", "Gilchrist_2021"))

#Rename columns
tabix_paths = dplyr::transmute(mappings, study, qtl_group,
                               tissue_ontology_id,
                               tissue_ontology_term,
                               tissue_label) %>%
  dplyr::left_join(condition_mapping, by = c("study", "qtl_group")) %>%
  dplyr::mutate(condition_label = ifelse(is.na(condition_label), "naive", condition_label)) %>%
  dplyr::mutate(tissue_ontology_term = paste0("\"", tissue_ontology_term, "\""))

#Process RNA-seq studies
rnaseq_paths = dplyr::anti_join(tabix_paths, microarray, by = "study") %>%
  tidyr::crossing(quants) %>%
  dplyr::mutate(ftp_path = paste("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats",
                                 study, quant_method,
                                 paste0(study, "_", quant_method, "_", qtl_group, ".all.tsv.gz"), sep ="/")) %>%
  dplyr::arrange(study, qtl_group, quant_method)
write.table(rnaseq_paths, "rnaseq_paths.tsv", sep = "\t", row.names = F, col.names = T, quote = F)


