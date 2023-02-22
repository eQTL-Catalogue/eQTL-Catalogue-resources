library("dplyr")

rnaseq_study_names <- c("Alasoo_2018", "BLUEPRINT_SE", "BLUEPRINT_PE", "BrainSeq", "GTEx", "FUSION", "GENCORD", "GEUVADIS", 
                        "HipSci", "Lepik_2017", "Nedelec_2016", "Quach_2016", "ROSMAP", "Schmiedel_2018", "Schwartzentruber_2018", 
                        "TwinsUK", "van_de_Bunt_2015", "Bossini-Castillo_2019", "CAP", "CommonMind", "Peng_2018", "PhLiPS", 
                        "iPSCORE","Braineac2","Steinberg_2020","Young_2019","Aygun_2021", "PISA","Walker_2019")
microarray_study_names <- c("CEDAR", "Fairfax_2014", "Kasela_2017", "Naranbhai_2015", "Fairfax_2012", "Gilchrist_2021")
protein_study_names <- c("Sun_2018")
file_names = c(rnaseq_study_names, microarray_study_names, protein_study_names)
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
message("Number of donors: ", length(unique(final_meta$genotype_id)))

#Calculate sample sizes
sample_sizes = dplyr::group_by(final_meta, study, qtl_group) %>% 
  dplyr::summarise(sample_size = length(sample_id)) %>% 
  dplyr::ungroup()

#Map to tissues
friendly_names = readr::read_tsv("ontology_mappings/friendly_names.tsv") %>%
  dplyr::select(tissue_ontology_id, tissue_label)
tissue_map = readr::read_tsv("ontology_mappings/tissue_ontology_mapping.tsv") %>%
  dplyr::select(study, qtl_group, tissue_ontology_id, tissue_ontology_term) %>%
  dplyr::left_join(friendly_names)
condition_labels = readr::read_tsv("ontology_mappings/condition_labels.tsv") 

dataset_table = dplyr::select(final_meta, study, study_file, qtl_group, protocol) %>% 
  dplyr::distinct() %>%
  dplyr::left_join(tissue_map) %>%
  dplyr::left_join(sample_sizes, by = c("study", "qtl_group")) %>%
  dplyr::left_join(condition_labels, by = c("study", "qtl_group")) %>%
  dplyr::mutate(condition_label = ifelse(is.na(condition_label), "naive", condition_label))
message("Number of cell types and tissues: ", length(unique(dataset_table$tissue_label)))

#Add quantification methods to datasets
quant_methods = dplyr::tibble(quant_method = c("ge","exon","tx","txrev","leafcutter"))
rnaseq_datasets = dplyr::filter(dataset_table, protocol %in% c("total", "poly(A)")) %>%
  dplyr::left_join(quant_methods, by = character())
microarray_datasets = dplyr::filter(dataset_table, protocol %in% c("HumanHT-12_V4")) %>%
  dplyr::mutate(quant_method = "microarray")
protein_datasets = dplyr::filter(dataset_table, protocol %in% c("SomaLogic")) %>%
  dplyr::mutate(quant_method = "aptamer")
all_datasets = dplyr::bind_rows(rnaseq_datasets, microarray_datasets, protein_datasets)

#Assing study and dataset ids
dataset_id_map = readr::read_tsv("data_tables/dataset_id_map.tsv")
dataset_metadata = dplyr::rename(all_datasets, study_label = study, 
                                 sample_group = qtl_group,
                                 tissue_id = tissue_ontology_id) %>% 
  dplyr::left_join(dataset_id_map, by = c("study_label", "quant_method", "sample_group")) %>%
  dplyr::select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size, quant_method) %>%
  dplyr::arrange(study_id, dataset_id)
write.table(dataset_metadata, "data_tables/dataset_metadata.tsv", sep = "\t", row.names = F, quote = F)


