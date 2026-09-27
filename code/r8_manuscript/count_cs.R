library("dplyr")

#Import R8 credbile sets
file_names = list.files("big_data/all_cs/", full.names = T)
ds_names = list.files("big_data/all_cs") %>% stringr::str_remove(".credible_sets.parquet")

file_list = setNames(as.list(file_names), ds_names)

cs_df = purrr::map_df(file_list, ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm, -rsid) %>%
                        dplyr::distinct() %>%
                        dplyr::group_by(cs_id) %>% dplyr::arrange(-pip) %>% 
                        dplyr::slice_head(n = 1), .id = "dataset_id")
cs_df = ungroup(cs_df)
arrow::write_parquet(cs_df, "big_data/r8_credible_set_leads.parquet")

ukkb_cs_df = arrow::read_parquet("big_data/UKBB_EUR_fine_mapping_with_meta_EUR_lead_variants.parquet")
ukbb_cs_size = dplyr::group_by(ukkb_cs_df, cs_id) %>% dplyr::summarise(cs_size = length(cs_id))

#Import r7 credible sets
file_names = list.files("big_data/r7_cs/", full.names = T)
ds_names = list.files("big_data/r7_cs") %>% stringr::str_remove(".credible_sets.tsv.gz")
file_list = setNames(as.list(file_names), ds_names)

cs_df = purrr::map_df(file_list, ~readr::read_tsv(.) %>%
                        dplyr::distinct() %>%
                        dplyr::group_by(cs_id) %>% dplyr::arrange(-pip) %>% 
                        dplyr::slice_head(n = 1), .id = "dataset_id")
cs_df = ungroup(cs_df)
r6_meta = readr::read_tsv("data_tables/dataset_metadata_r6.tsv")
r6_cs_df = dplyr::filter(cs_df, dataset_id %in% r6_meta$dataset_id)
arrow::write_parquet(r6_cs_df, "big_data/r6_credible_set_leads.parquet")


#Compare credible set count for r6 and r8
r6_cs_df = arrow::read_parquet("big_data/r6_credible_set_leads.parquet")
r8_cs_df = arrow::read_parquet("big_data/r8_credible_set_leads.parquet")

#How many credible sets originate from datasets not present in r6?
dplyr::anti_join(r8_cs_df, r6_meta, by = "dataset_id")

#How many credible sets come from single-cell datasets
r8_meta = readr::read_tsv("data_tables/dataset_metadata_r8.tsv")
r8_sceqtl = dplyr::filter(r8_meta, study_type == "single-cell")
dplyr::semi_join(r8_cs_df, r8_sceqtl, by = "dataset_id")

#How many are from IBDverse?
r8_ibd = dplyr::filter(r8_meta, study_label == "IBDverse")
dplyr::semi_join(r8_cs_df, r8_ibd, by = "dataset_id")

#How many are majiq credible sets?
r8_majiq = dplyr::filter(r8_meta, quant_method == "majiq")
dplyr::semi_join(r8_cs_df, r8_majiq, by = "dataset_id")

#How many novel credible sets do we detect for GTEx?
r6_gtex = dplyr::filter(r6_meta, study_label == "GTEx")
dplyr::semi_join(r6_cs_df, r6_gtex, by = "dataset_id")

r8_gtex = dplyr::filter(r8_meta, study_label == "GTEx_v10", quant_method != "majiq")
dplyr::semi_join(r8_cs_df, r8_gtex, by = "dataset_id")


#How many novel eqtl credible sets do we detect for GTEx?
r6_gtex = dplyr::filter(r6_meta, study_label == "GTEx", quant_method == "ge")
dplyr::semi_join(r6_cs_df, r6_gtex, by = "dataset_id")

r8_gtex = dplyr::filter(r8_meta, study_label == "GTEx_v10", quant_method == "ge")
dplyr::semi_join(r8_cs_df, r8_gtex, by = "dataset_id")

#How many non-GTEx novel credible sets do we detect in r8
r6_not_gtex = dplyr::filter(r6_meta, study_label != c("GTEx"))
r6_not = dplyr::semi_join(r6_cs_df, r6_not_gtex, by = "dataset_id")
r6_cs_count = dplyr::group_by(r6_not, dataset_id) %>% dplyr::summarise(r6_cs_count = n())

r8_not = dplyr::semi_join(r8_cs_df, r6_not_gtex, by = "dataset_id")
r8_cs_count = dplyr::group_by(r8_not, dataset_id) %>% dplyr::summarise(r8_cs_count = n())

percent_increase = dplyr::left_join(r6_cs_count, r8_cs_count) %>% 
  dplyr::left_join(r6_meta) %>%
  dplyr::mutate(pct_increase = (r8_cs_count - r6_cs_count)/r6_cs_count)
ge_increase = dplyr::filter(percent_increase, quant_method %in% c("ge", "microarray"))
exon_increase = dplyr::filter(percent_increase, quant_method == "exon")
other_increase = dplyr::filter(percent_increase, !(quant_method %in% c("ge", "microarray", "exon")))


#Lepik_2017
lepik_cs_df = arrow::read_parquet("big_data/all_cs/QTD000373.credible_sets.parquet")
mage_cs_df = arrow::read_parquet("big_data/all_cs/QTD000764.credible_sets.parquet")
geuvadis_cs_df = arrow::read_parquet("big_data/all_cs/QTD000110.credible_sets.parquet")


#Infer samples sizes for GTEx credible sets
gtex_metadata = readr::read_tsv("data_tables/dataset_metadata_r8.tsv") %>%
  dplyr::filter(study_label == "GTEx_v10", quant_method == "ge")

gtex_file_list = file_list[gtex_metadata$dataset_id]

#Import gtex credible sets
gtex_cs = purrr::map_df(gtex_file_list, ~arrow::read_parquet(.) %>% dplyr::select(-r2, -median_tpm) %>%
                        dplyr::distinct(), .id = "dataset_id")
gtex_v10_sample_sizes = dplyr::mutate(cs_by_ds, inferred_n = round(ac/maf)/2) %>% 
  dplyr::select(dataset_id, inferred_n) %>% 
  dplyr::group_by(dataset_id) %>% 
  dplyr::summarise(sample_size = get_mode(inferred_n)) %>% dplyr::arrange(dataset_id)
write.table(gtex_v10_sample_sizes, "data_tables/gtex_v10_sample_size.tsv", sep = "\t", row.names = F, quote = F)


