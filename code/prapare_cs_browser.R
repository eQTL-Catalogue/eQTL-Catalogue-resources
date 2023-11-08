library("dplyr")
library("arrow")
library("stringr")

#Import dataset metadata
dataset_metadata = readr::read_tsv("../eQTL-Catalogue-resources/data_tables/dataset_metadata.tsv") %>%
  dplyr::mutate(dataset_label = paste(study_label, tissue_label, condition_label, sep = ":")) 
dataset_labels = dplyr::select(dataset_metadata, dataset_id, dataset_label, quant_method)

#Make a list of paths
paths = file.path("large_data/credible_sets_r6/", paste0(dataset_metadata$dataset_id, ".credible_sets.tsv.gz"))
path_list = setNames(paths, dataset_metadata$dataset_id)

#Import tables
cs_table = purrr::map_df(path_list, readr::read_tsv, .id = "dataset_id")
cs_table = readRDS("large_data/credible_sets_r6.rds")

#Load gene name mapping
gene_meta = readr::read_tsv("https://zenodo.org/record/7808390/files/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz") %>%
  dplyr::select(gene_id, gene_name)
microarray_meta = readr::read_tsv("https://zenodo.org/record/7808390/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::anti_join(gene_meta, by = "gene_id")
somalogic_meta = readr::read_tsv("https://zenodo.org/record/7808390/files/SomaLogic_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::anti_join(gene_meta, by = "gene_id")
gene_names = dplyr::bind_rows(gene_meta, microarray_meta, somalogic_meta) %>% dplyr::distinct()

#Assign unique credible set ids
cs_table_new = dplyr::mutate(cs_table, credible_set = paste(dataset_id, cs_id, sep = "_")) %>%
  dplyr::left_join(dataset_labels, by = "dataset_id") %>%
  dplyr::left_join(gene_names, by = "gene_id") %>%
  dplyr::transmute(molecular_trait_id, gene_id, gene_name, credible_set,
                   variant, rs_id = rsid, credible_set_size = cs_size, pip, p_value = pvalue, beta, se, dataset_id, dataset_label, quant_method, z, region)
arrow::write_parquet(cs_table_new, "large_data/eQTL_Catalogue_r6_credible_sets.parquet")

#Add credible set ids to coverage plot metadata
cs_ids = dplyr::select(cs_table_new, molecular_trait_id, variant, dataset_id, credible_set) %>%
  dplyr::distinct()
plot_meta = arrow::read_parquet("large_data/coverage_plot_metadata.parquet") %>% dplyr::rename(dataset_id = dataset)

#Rename leafcutter molecular_trait_ids
lc_datasets = dplyr::filter(dataset_metadata, quant_method == "leafcutter")
lc_plots = dplyr::semi_join(plot_meta, lc_datasets, by = "dataset_id") %>%
  dplyr::mutate(molecular_trait_id = str_replace(molecular_trait_id, "_", ":") %>% str_replace("_", ":") %>% str_replace("_", ":")) %>%
  dplyr::left_join(cs_ids, by = c("molecular_trait_id", "variant", "dataset_id"))
non_lc_plots = dplyr::anti_join(plot_meta, lc_datasets, by = "dataset_id") %>%
  dplyr::left_join(cs_ids, by = c("molecular_trait_id", "variant", "dataset_id"))
all_plots = dplyr::bind_rows(lc_plots, non_lc_plots)

#Write new plot metadata to disk
arrow::write_parquet(all_plots, "large_data/eQTL_Catalogue_r6_coverage_plots.parquet")




tabix_table = readr::read_tsv("../tabix/tabix_ftp_paths.tsv") %>%
  dplyr::filter(quant_method %in% c("ge","microarray")) %>%
  dplyr::mutate(study_proxy = ifelse(study == "BLUEPRINT", "BLUEPRINT_SE", study)) %>%
  dplyr::mutate(study_proxy = ifelse(study_proxy == "BLUEPRINT_SE" & qtl_group == "T-cell", "BLUEPRINT_PE", study_proxy)) %>%
  dplyr::mutate(cs_path = paste0(study_proxy, "_",quant_method,"_", qtl_group, ".credible_sets.tsv.gz")) %>%
  dplyr::mutate(dataset = paste(study, tissue_label, condition_label, sep = ":"))

#Make a list of paths
paths = file.path("../data/susie_merged", tabix_table$cs_path)
path_list = setNames(paths, tabix_table$dataset)


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

