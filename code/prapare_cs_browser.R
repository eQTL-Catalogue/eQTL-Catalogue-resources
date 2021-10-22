library("dplyr")

#Import dataset metadata
tabix_table = readr::read_tsv("../tabix/tabix_ftp_paths.tsv") %>%
  dplyr::filter(quant_method %in% c("ge","microarray")) %>%
  dplyr::mutate(study_proxy = ifelse(study == "BLUEPRINT", "BLUEPRINT_SE", study)) %>%
  dplyr::mutate(study_proxy = ifelse(study_proxy == "BLUEPRINT_SE" & qtl_group == "T-cell", "BLUEPRINT_PE", study_proxy)) %>%
  dplyr::mutate(cs_path = paste0(study_proxy, "_",quant_method,"_", qtl_group, ".credible_sets.tsv.gz")) %>%
  dplyr::mutate(dataset = paste(study, tissue_label, condition_label, sep = ":"))

#Make a list of paths
paths = file.path("../data/susie_merged", tabix_table$cs_path)
path_list = setNames(paths, tabix_table$dataset)

#Import tables
dataset_table = purrr::map_df(path_list, readr::read_tsv, .id = "dataset")

#Calculate dataset indices
ds_indices = dplyr::mutate(tabix_table, dataset_index = paste0("D", c(1:length(dataset)))) %>%
  dplyr::select(dataset, dataset_index)

#Import gene and probe metadata
#Import gene metadata
gene_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::select(phenotype_id, gene_name)
probe_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::select(phenotype_id, gene_name)
gene_map = dplyr::bind_rows(gene_meta, probe_meta) %>%
  dplyr::rename(molecular_trait_id = phenotype_id)

#Make the final table
final_table = dplyr::left_join(dataset_table, ds_indices, by = "dataset") %>%
  dplyr::left_join(gene_map, by = "molecular_trait_id") 

save_table = dplyr::mutate(final_table, credible_set = paste0(cs_id, "_", dataset_index)) %>%
  dplyr::rename(credible_set_size = cs_size) %>%
  dplyr::select(molecular_trait_id, gene_name, credible_set, variant, rsid, credible_set_size, pip, pvalue, beta, se, dataset) %>%
  dplyr::rename(rs_id = rsid, p_value = pvalue)
file_handle = gzfile("../data/credible_set_table.tsv.gz","w")
write.table(save_table, file_handle, sep = "\t", row.names = F, col.names = T, quote = FALSE)
close(file_handle)
