library("dplyr")

#Import data
r7 = readr::read_tsv("data_tables/dataset_metadata_r7.tsv")

#Merge credible sets by quantification method
file_names = list.files("big_data/r7_cs/", full.names = T)
ds_names = list.files("big_data/r7_cs") %>% stringr::str_remove(".credible_sets.tsv.gz")
file_list = setNames(as.list(file_names), ds_names)

#Remove missing cs file form r7 metadata
r7 = dplyr::filter(r7, dataset_id != "QTD000736")

#Import ge/microarray/aptamer credible sets
ge_cs_meta = dplyr::filter(r7, quant_method %in% c("ge","microarray","aptamer"))
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~readr::read_tsv(.) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = tidyr::separate(ge_cs, variant, c("chromosome", "position", "ref", "alt"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(chromosome = stringr::str_replace(chromosome, "chr", "")) %>%
  dplyr::mutate(position = as.integer(position)) %>%
  dplyr::arrange(chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant_r7/eQTL_Catalogue_r7-beta_cs_ge_071026.parquet")

#Import leafcutter credible sets
ge_cs_meta = dplyr::filter(r7, quant_method %in% c("leafcutter"))
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~readr::read_tsv(.) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = tidyr::separate(ge_cs, variant, c("chromosome", "position", "ref", "alt"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(chromosome = stringr::str_replace(chromosome, "chr", "")) %>%
  dplyr::mutate(position = as.integer(position)) %>%
  dplyr::arrange(chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant_r7/eQTL_Catalogue_r7_cs_leafcutter_071026.parquet")

#Import exon credible sets
ge_cs_meta = dplyr::filter(r7, quant_method %in% c("exon"))
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~readr::read_tsv(.) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = tidyr::separate(ge_cs, variant, c("chromosome", "position", "ref", "alt"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(chromosome = stringr::str_replace(chromosome, "chr", "")) %>%
  dplyr::mutate(position = as.integer(position)) %>%
  dplyr::arrange(chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant_r7/eQTL_Catalogue_r7_cs_exon_071026.parquet")

#Import tx credible sets
ge_cs_meta = dplyr::filter(r7, quant_method %in% c("tx"))
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~readr::read_tsv(.) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = tidyr::separate(ge_cs, variant, c("chromosome", "position", "ref", "alt"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(chromosome = stringr::str_replace(chromosome, "chr", "")) %>%
  dplyr::mutate(position = as.integer(position)) %>%
  dplyr::arrange(chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant_r7/eQTL_Catalogue_r7_cs_tx_071026.parquet")

#Import txrev credible sets
ge_cs_meta = dplyr::filter(r7, quant_method %in% c("txrev"))
ge_cs = purrr::map_df(file_list[ge_cs_meta$dataset_id], ~readr::read_tsv(.) %>%
                        dplyr::distinct(), .id = "dataset_id")

ge_cs_sorted = tidyr::separate(ge_cs, variant, c("chromosome", "position", "ref", "alt"), sep = "_", remove = FALSE) %>%
  dplyr::mutate(chromosome = stringr::str_replace(chromosome, "chr", "")) %>%
  dplyr::mutate(position = as.integer(position)) %>%
  dplyr::arrange(chromosome, position, dataset_id, molecular_trait_id) %>%
  dplyr::mutate(cs_id = paste0(dataset_id, "_", cs_id))
arrow::write_parquet(ge_cs_sorted, "big_data/cs_by_quant_r7/eQTL_Catalogue_r7_cs_txrev_071026.parquet")

