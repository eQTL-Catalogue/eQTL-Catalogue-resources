library("dplyr")

#Import dataset metadata
dataset_metadata = readr::read_tsv("data_tables/dataset_metadata.tsv")

#Add ftp paths
tabix_data = dplyr::mutate(dataset_metadata, 
                           ftp_path = paste0("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/", study_id, "/", dataset_id, "/", dataset_id)) %>%
  dplyr::mutate(ftp_path = ifelse(quant_method %in% c("ge", "microarray", "aptamer"), paste0(ftp_path, ".all.tsv.gz"),paste0(ftp_path, ".cc.tsv.gz"))) %>%
  dplyr::mutate(ftp_cs_path = paste0("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/", study_id, "/", dataset_id, "/", dataset_id, ".credible_sets.tsv.gz")) %>%
  dplyr::mutate(ftp_lbf_path = paste0("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/", study_id, "/", dataset_id, "/", dataset_id, ".lbf_variable.txt.gz"))
write.table(tabix_data, "tabix/tabix_ftp_paths.tsv", sep = "\t", row.names = F, col.names = T, quote = F)


