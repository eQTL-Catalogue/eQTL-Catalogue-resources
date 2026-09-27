library("dplyr")
r7_meta = readr::read_tsv("data_tables/dataset_metadata_r7.tsv")
r6_meta = dplyr::filter(r7_meta, !(study_id %in% c("QTS000032", "QTS000033", "QTS000034", "QTS000042", 
                                                   "QTS000036","QTS000037", "QTS000036",
                                                   "QTS000037", "QTS000038", "QTS000039", "QTS000040",
                                                   "QTS000041")))
write.table(r6_meta, "data_tables/dataset_metadata_r6.tsv", quote = F, row.names = F, sep = "\t")
r8_meta = readr::read_tsv("data_tables/dataset_metadata_r8.tsv")

#Count unique bulk studies
dplyr::select(r6_meta, study_id) %>% distinct() #32
dplyr::filter(r8_meta, study_type == "bulk") %>% dplyr::select(study_id) %>% distinct() #43

#Count unique bulk datasets
dplyr::select(r6_meta, dataset_id) %>% distinct() #564
dplyr::filter(r8_meta, study_type == "bulk") %>% dplyr::select(dataset_id) %>% distinct() #884

#Count sample groups
r6_meta %>% dplyr::select(study_id, sample_group) %>% dplyr::distinct()
dplyr::filter(r8_meta, study_type == "bulk") %>% dplyr::select(study_id, sample_group) %>% dplyr::distinct()

#Compare GTEx v8 (r6) and GTEx v10 (r8) sample size
v8_sample_size = dplyr::filter(r6_meta, study_id == "QTS000015", quant_method == "ge") %>% 
  dplyr::transmute(dataset_id, sample_group, v8_sample_size = sample_size)
v10_sample_size = dplyr::filter(r8_meta, study_id == "QTS000015", quant_method == "ge") %>% 
  transmute(dataset_id, v10_sample_size = sample_size)
gtex_v10 = dplyr::left_join(v8_sample_size, v10_sample_size) %>%
  dplyr::mutate(percent_increase = (v10_sample_size-v8_sample_size)/v8_sample_size)
write.table(gtex_v10, "manuscript/r8_tables/GTEx_v8_vs_v10_sample_sizes.tsv", sep = "\t", quote = F, row.names = F)


