library("dplyr")

#Import sceQTL sample sizes
sample_sizes = readr::read_table("data_tables/sceQTL_num_individuals_per_dataset.tsv")

#Import dataset id map
dataset_id_map = readr::read_tsv("data_tables/dataset_id_map.tsv")

#Map to tissues
friendly_names = readr::read_tsv("ontology_mappings/friendly_names.tsv") %>%
  dplyr::select(tissue_ontology_id, tissue_label)
tissue_map = readr::read_tsv("ontology_mappings/tissue_ontology_mapping.tsv") %>%
  dplyr::select(study, qtl_group, tissue_ontology_id, tissue_ontology_term) %>%
  dplyr::left_join(friendly_names) %>%
  dplyr::rename(sample_group = qtl_group) %>%
  dplyr::rename(study_label = study) %>%
  dplyr::rename(tissue_id = tissue_ontology_id)
condition_labels = readr::read_tsv("ontology_mappings/condition_labels.tsv") %>%
  dplyr::rename(sample_group = qtl_group) %>%
  dplyr::rename(study_label = study)

dataset_table = dplyr::left_join(sample_sizes, dataset_id_map, by = "dataset_id") %>%
  dplyr::left_join(tissue_map, by = c("study_label", "sample_group")) %>%
  dplyr::left_join(condition_labels, by = c("study_label", "sample_group")) %>%
  dplyr::mutate(condition_label = ifelse(is.na(condition_label), "naive", condition_label)) %>% 
  dplyr::select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size, quant_method) %>%
  dplyr::arrange(study_id, dataset_id)

write.table(dataset_table, "data_tables/dataset_metadata_sceQTL.tsv", sep = "\t", row.names = F, quote = F)