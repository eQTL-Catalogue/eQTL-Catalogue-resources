library("dplyr")

rnaseq_study_names <- c("Alasoo_2018", "BLUEPRINT_SE", "BLUEPRINT_PE", "BrainSeq", "GTEx", "FUSION", "GENCORD", "GEUVADIS", 
                        "HipSci", "Lepik_2017", "Nedelec_2016", "Quach_2016", "ROSMAP", "Schmiedel_2018", "Schwartzentruber_2018", 
                        "TwinsUK", "van_de_Bunt_2015", "Bossini-Castillo_2019", "CAP", "CommonMind", "Peng_2018", "PhLiPS", 
                        "iPSCORE","Braineac2","Steinberg_2020","Young_2019")
microarray_study_names <- c("CEDAR", "Fairfax_2014", "Kasela_2017", "Naranbhai_2015", "Fairfax_2012", "Gilchrist_2021")
file_names = c(rnaseq_study_names, microarray_study_names)
meta_path = "../SampleArcheology/studies/cleaned/"


mandatory_columns_meta <- c("sample_id", "genotype_id", "sex", "cell_type","condition", "qtl_group", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")

#Import all studies
merged_meta <- dplyr::tibble()
for (study in file_names) {
  study_meta <- readr::read_tsv(stringr::str_interp("${meta_path}/${study}.tsv"))
  study_meta_filt <- study_meta[,mandatory_columns_meta]
  merged_meta <- merged_meta %>% rbind(study_meta_filt)
}

#Fix BLUEPRINT study name
final_meta = dplyr::mutate(merged_meta, study_file = study) %>%
  dplyr::mutate(study = ifelse(study %in% c("BLUEPRINT_PE", "BLUEPRINT_SE"), "BLUEPRINT", study)) %>%
  dplyr::filter(rna_qc_passed, genotype_qc_passed) %>%
  dplyr::mutate(dataset_name = paste0(study_file, "_", qtl_group))

message("Number of studies: ", length(unique(final_meta$study)))
message("Number of datasets: ", length(unique(final_meta$dataset_name)))
message("Number of donors:", length(unique(final_meta$genotype_id)))

#Map to tissues
friendly_names = readr::read_tsv("ontology_mappings/friendly_names.tsv") %>%
  dplyr::select(tissue_ontology_id, tissue_label)
tissue_map = readr::read_tsv("ontology_mappings/tissue_ontology_mapping.tsv") %>%
  dplyr::select(study, qtl_group, tissue_ontology_id, tissue_ontology_term) %>%
  dplyr::left_join(friendly_names)

dataset_table = dplyr::select(final_meta, study, study_file, qtl_group) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(tissue_map)
message("Number of cell types and tissues: ", length(unique(dataset_table$tissue_label)))

#Assign study ids
studies = dplyr::select(dataset_table, study) %>% 
  dplyr::distinct() %>% dplyr::arrange(study)
study_ids = dplyr::mutate(studies, study_index = c(1:length(study))) %>% 
  dplyr::mutate(study_id = ifelse(study_index < 10, paste0("QTS", "00000", study_index), paste0("QTS", "0000", study_index))) %>%
  dplyr::select(study, study_id)
write.table(study_ids, "data_tables/study_id_map.tsv", sep = "\t", row.names = F, quote = F)

#Assign dataset ids
dataset_ids = dplyr::left_join(dataset_table, study_ids) %>%
  dplyr::arrange(study) %>%
  dplyr::mutate(study_index = c(1:length(study))) %>%
  dplyr::mutate(padded_study_idx = stringr::str_pad(study_index, 6, "0", side = "left")) %>%
  dplyr::mutate(dataset_id = paste0("QTD", padded_study_idx)) %>%
  dplyr::select(study, qtl_group, study_id, dataset_id)
write.table(dataset_ids, "data_tables/dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)
