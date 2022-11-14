
makeDatasetMetadata <- function(study_label, sample_groups, quant_methods = c("ge","exon","tx","txrev","leafcutter"), study_index_start, dataset_index_start){
  
  quant_methods = dplyr::tibble(quant_method = quant_methods)
  
  df = dplyr::tibble(study_label = study_label, sample_group = sample_groups) %>%
    dplyr::left_join(quant_methods, by = character()) %>%
    dplyr::mutate(study_index = study_index_start) %>%
    dplyr::mutate(dataset_index = c(dataset_index_start:(dataset_index_start+length(study_label)-1))) %>%
    dplyr::mutate(padded_study_idx = stringr::str_pad(study_index, 6, "0", side = "left")) %>%
    dplyr::mutate(padded_dataset_idx = stringr::str_pad(dataset_index, 6, "0", side = "left")) %>%
    dplyr::mutate(study_id = paste0("QTS", padded_study_idx)) %>%
    dplyr::mutate(dataset_id = paste0("QTD", padded_dataset_idx)) %>%
    dplyr::select(study_label,sample_group,quant_method,study_id,dataset_id)
  return(df)
}

id_map = makeDatasetMetadata("Aygun_2021", c("Progenitor", "Neuron"), study_index_start = 32, dataset_index_start = 564)
write.table(id_map, "data_tables/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)
