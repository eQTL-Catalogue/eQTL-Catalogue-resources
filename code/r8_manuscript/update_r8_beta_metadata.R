library("dplyr")

#Import data
r8 = readr::read_tsv("data_tables/dataset_metadata_r8.tsv")
r8_beta = readr::read_tsv("data_tables/dataset_metadata_r8_beta.tsv")

#Find newly added datasets
added = dplyr::filter(r8, study_label %in% c("GTEx_v10", "MAGE", "IBDverse"))

r8_beta_new = dplyr::bind_rows(r8_beta, added) %>%
  dplyr::arrange(study_id, dataset_id)
write.table(r8_beta_new, "data_tables/dataset_metadata_r8_beta.tsv", quote = F, row.names = F, sep = "\t")

#Update GTEx_v10 sample_size
gtex_v10 = dplyr::filter(r8, study_label == "GTEx_v10")

gtex_v10_ss = readr::read_tsv("data_tables/gtex_v10_sample_size.tsv") %>% 
  dplyr::left_join(dplyr::select(gtex_v10, dataset_id, sample_group)) %>%
  dplyr::select(-dataset_id)
gtex_v10_updated_ss = dplyr::select(gtex_v10, -sample_size) %>% 
  dplyr::left_join(gtex_v10_ss) %>% 
  dplyr::select(study_id,dataset_id,study_label,sample_group, tissue_id, tissue_label, condition_label, sample_size, everything())
write.table(gtex_v10_updated_ss, "data_tables/gtex_v10_update_sample_sizes.tsv", sep = "\t", quote = F, row.names = F)

#Merge credible sets by quantification method
file_names = list.files("big_data/all_cs/", full.names = T)
ds_names = list.files("big_data/all_cs") %>% stringr::str_remove(".credible_sets.parquet")
file_list = setNames(as.list(file_names), ds_names)

#Import ge credible sets
ge_cs_meta = dplyr::filter(r8_beta_new, quant_method == "ge")
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                                dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = dplyr::arrange(ge_cs, chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant/eQTL_Catalogue_r8-beta_cs_ge_190926.parquet")

#Leafcutter
ge_cs_meta = dplyr::filter(r8_beta_new, quant_method == "leafcutter")
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = dplyr::arrange(ge_cs, chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant/eQTL_Catalogue_r8-beta_cs_leafcutter_190926.parquet")

#MAJIQ
ge_cs_meta = dplyr::filter(r8_beta_new, quant_method == "majiq")
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = dplyr::arrange(ge_cs, chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant/eQTL_Catalogue_r8-beta_cs_majiq_190926.parquet")

#Exon
ge_cs_meta = dplyr::filter(r8_beta_new, quant_method == "exon")
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = dplyr::arrange(ge_cs, chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant/eQTL_Catalogue_r8-beta_cs_exon_190926.parquet")

#Tx
ge_cs_meta = dplyr::filter(r8_beta_new, quant_method == "tx")
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = dplyr::arrange(ge_cs, chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant/eQTL_Catalogue_r8-beta_cs_tx_190926.parquet")

#Txrev
ge_cs_meta = dplyr::filter(r8_beta_new, quant_method == "txrev")
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = dplyr::arrange(ge_cs, chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant/eQTL_Catalogue_r8-beta_cs_txrev_190926.parquet")


