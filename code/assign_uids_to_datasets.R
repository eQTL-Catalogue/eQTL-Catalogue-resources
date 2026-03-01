library("dplyr")
makeDatasetMetadata <- function(study_label, sample_groups, quant_methods = c("ge","exon","tx","txrev","leafcutter", "majiq"), study_index_start, dataset_index_start){
  
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

aygun_map = makeDatasetMetadata("Aygun_2021", c("Progenitor", "Neuron"), study_index_start = 32, dataset_index_start = 564)
pisa_map = makeDatasetMetadata("PISA", c("pancreatic_islet"), study_index_start = 33, dataset_index_start = 574)
walker_map = makeDatasetMetadata("Walker_2019", c("Neocortex"), study_index_start = 34, dataset_index_start = 579)
sun_map = makeDatasetMetadata("Sun_2018", c("plasma"), quant_methods = "somalogic", study_index_start = 35, dataset_index_start = 584)
new_meta = dplyr::bind_rows(aygun_map, pisa_map, walker_map, sun_map)

write.table(new_meta, "data_tables/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)

#Add single-cell datasets
randolph_map = makeDatasetMetadata("Randolph_2021", c("B_NI","B_flu","CD4_T_NI","CD4_T_flu","CD8_T_NI","CD8_T_flu","NK_NI","NK_flu","highly_infected_flu","infected_monocytes_flu","monocytes_NI","monocytes_flu"), quant_methods = "ge", study_index_start = 36, dataset_index_start = 585)
perez_map = makeDatasetMetadata("Perez_2022", c("B","NK","Prolif","T4","T8","cDC","cM","ncM","pDC"), quant_methods = "ge", study_index_start = 37, dataset_index_start = 597)
onek1k_map = makeDatasetMetadata("OneK1K", c("B_intermediate","B_memory","B_naive","CD14_Mono","CD16_Mono","CD4_CTL","CD4_Naive","CD4_TCM","CD4_TEM","CD8_Naive","CD8_TCM","CD8_TEM","HSPC","MAIT","NK","NK_CD56bright","NK_Proliferating","Plasmablast","Platelet","Treg","cDC2","dnT","gdT","pDC"), quant_methods = "ge", study_index_start = 38, dataset_index_start = 606)
new_meta = dplyr::bind_rows(randolph_map, perez_map, onek1k_map)
write.table(new_meta, "data_tables/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)

#Add Jerber_2021, Nathan_2022, Cytoimmgen
jerber_celltypes = readr::read_tsv("~/Downloads/celltype_cell_counts.Jerber_2021.tsv")
jerber_map = makeDatasetMetadata("Jerber_2021", sample_groups = jerber_celltypes$celltype, quant_methods = "ge", study_index_start = 39, dataset_index_start = 630)
nathan_celltypes = readr::read_tsv("~/Downloads/celltype_cell_counts.Nathan_2022.tsv")
nathan_map = makeDatasetMetadata("Nathan_2022", sample_groups = nathan_celltypes$celltype, quant_methods = "ge", study_index_start = 40, dataset_index_start = 660)
cytoimmgen_celltypes = readr::read_tsv("~/Downloads/celltype_cell_counts.Cytoimmgen.tsv")
cytoimmgen_map = makeDatasetMetadata("Cytoimmgen", sample_groups = cytoimmgen_celltypes$celltype, quant_methods = "ge", study_index_start = 41, dataset_index_start = 689)
new_meta = dplyr::bind_rows(jerber_map, nathan_map, cytoimmgen_map)
write.table(new_meta, "data_tables/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)

#Add MAJIQ ids for GTEx
ds_meta = readr::read_tsv("~/projects/eQTL-Catalogue-resources/data_tables/dataset_metadata_upcoming.tsv")
gtex_tissues = dplyr::filter(ds_meta, study_label == "GTEx", quant_method == "ge")

gtex_majiq_dsids = makeDatasetMetadata("GTEx", gtex_tissues$sample_group, quant_methods = c("majiq"), 15, 951)
write.table(gtex_majiq_dsids, "~/Downloads/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)

tissue_meta = dplyr::select(gtex_tissues, sample_group, tissue_id, tissue_label, condition_label, sample_size, pmid, study_type)

gtex_majiq_meta = dplyr::left_join(gtex_majiq_dsids, tissue_meta, by = "sample_group") %>%
  dplyr::select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size, quant_method, pmid, study_type)
write.table(gtex_majiq_meta, "~/Downloads/gtex_meta.tsv", sep = "\t", row.names = F, quote = F)

#Generate ids for INTERVAL RNA
interval_rna = makeDatasetMetadata("INTERVAL_RNA", "blood", quant_methods = c("ge","leafcutter","majiq"), 53, 1000)
interval_rna_wgs = makeDatasetMetadata("INTERVAL_RNA_WGS", "blood", quant_methods = c("ge","leafcutter","majiq"), 54, 1003)
interval_meta = dplyr::bind_rows(interval_rna, interval_rna_wgs) %>%
  dplyr::transmute(study_id, dataset_id, study_label, sample_group, tissue_id = "UBERON_0000178", tissue_label = "blood",
                   condition_label = "naive", sample_size = NA, quant_method, pmid = "40038547", study_type = "bulk")
write.table(interval_meta, "~/Downloads/interval_meta.tsv", sep = "\t", row.names = F, quote = F)

#AFR_LCL
afr_lcl = makeDatasetMetadata("AFR_LCL", "LCL", study_index_start = 55, dataset_index_start = 1006)
afr_lcl_meta = afr_lcl %>%
  dplyr::transmute(study_id, dataset_id, study_label, sample_group, tissue_id = "EFO_0005292", tissue_label = "LCL",
                   condition_label = "naive", sample_size = NA, quant_method, pmid = "NA", study_type = "bulk")
write.table(afr_lcl_meta, "~/Downloads/afr_lcl_meta.tsv", sep = "\t", row.names = F, quote = F)

#Alasoo_2018 majiq
ds_meta = readr::read_tsv("~/projects/eQTL-Catalogue-resources/data_tables/dataset_metadata_r7.tsv")
Alasoo_2018_tissues = dplyr::filter(ds_meta, study_label == "Alasoo_2018", quant_method == "ge")

alasoo_2018 = makeDatasetMetadata("Alasoo_2018", Alasoo_2018_tissues$sample_group, quant_methods = c("majiq"), 1, 1026)
write.table(alasoo_2018, "~/Downloads/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)

tissue_meta = dplyr::select(Alasoo_2018_tissues, sample_group, tissue_id, tissue_label, condition_label, sample_size, pmid, study_type)

alasoo_majiq_meta = dplyr::left_join(alasoo_2018, tissue_meta, by = "sample_group") %>%
  dplyr::select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size, quant_method, pmid, study_type)
write.table(alasoo_majiq_meta, "~/Downloads/gtex_meta.tsv", sep = "\t", row.names = F, quote = F)

#Nassiri_2025
nassiri_sample_size = readr::read_tsv("../../SampleArcheology/studies/cleaned/Nassiri_2025.tsv") %>%
  dplyr::filter(genotype_qc_passed, rna_qc_passed) %>%
  dplyr::group_by(qtl_group) %>%
  dplyr::summarise(sample_size = n())
nassiri_conditions = readr::read_tsv("../ontology_mappings/condition_labels.tsv")
friendly_names = readr::read_tsv("../ontology_mappings/friendly_names.tsv")

nassiri_mappings = readr::read_tsv("../ontology_mappings/tissue_ontology_mapping.tsv") %>% 
  dplyr::filter(study == "Nassiri_2025") %>%
  dplyr::left_join(nassiri_sample_size) %>%
  dplyr::left_join(nassiri_conditions) %>%
  dplyr::left_join(friendly_names) %>%
  dplyr::rename(sample_group = qtl_group) %>%
  dplyr::rename(study_label = study) %>%
  dplyr::rename(tissue_id = tissue_ontology_id)

nassiri_dataset_ids = makeDatasetMetadata("Nassiri_2025", nassiri_mappings$sample_group, quant_methods = c("ge","exon","tx","txrev","leafcutter","majiq"), 58, 1030)
nassiri_meta = dplyr::left_join(nassiri_dataset_ids, nassiri_mappings, by = c("study_label", "sample_group")) %>%
  dplyr::transmute(study_id, dataset_id, study_label, sample_group, tissue_id, 
                   tissue_label, condition_label, sample_size, quant_method, pmid = "41022800", study_type = "bulk")
write.table(nassiri_meta, "~/Downloads/Nassiri_meta.tsv", sep = "\t", row.names = F, quote = F)


#### MacroMap ####
#Import sample metadata
macromap_meta = readr::read_tsv("../SampleArcheology/studies/cleaned/MacroMap.tsv")
macromap_sample_size = macromap_meta %>%
  dplyr::filter(genotype_qc_passed, rna_qc_passed) %>%
  dplyr::group_by(qtl_group) %>%
  dplyr::summarise(sample_size = n())
conditions = readr::read_tsv("ontology_mappings/condition_labels.tsv")
friendly_names = readr::read_tsv("ontology_mappings/friendly_names.tsv")

macropmap_mappings = readr::read_tsv("ontology_mappings/tissue_ontology_mapping.tsv") %>% 
  dplyr::filter(study == "MacroMap") %>%
  dplyr::left_join(macromap_sample_size) %>%
  dplyr::left_join(conditions) %>%
  dplyr::left_join(friendly_names) %>%
  dplyr::rename(sample_group = qtl_group) %>%
  dplyr::rename(study_label = study) %>%
  dplyr::rename(tissue_id = tissue_ontology_id)

macromap_dataset_ids = makeDatasetMetadata("MacroMap", macromap_sample_size$qtl_group, quant_methods = c("ge","exon","tx","txrev","leafcutter","majiq"), 59, 1055)
write.table(macromap_dataset_ids, "~/Downloads/MacroMap_ids.tsv", sep = "\t", row.names = F, quote = F)

macromap_ds_meta = dplyr::left_join(macromap_dataset_ids, macropmap_mappings, by = c("study_label", "sample_group")) %>%
  dplyr::transmute(study_id, dataset_id, study_label, sample_group, tissue_id, 
                   tissue_label, condition_label, sample_size, quant_method, pmid = "40866338", study_type = "bulk")
write.table(macromap_ds_meta, "~/Downloads/MacroMap_meta.tsv", sep = "\t", row.names = F, quote = F)


#Add MAJIQ ids for Schmiedel_2018
#Add MAJIQ ids for GTEx
ds_meta = readr::read_tsv("~/projects/eQTL-Catalogue-resources/data_tables/dataset_metadata_r7.tsv")
gtex_tissues = dplyr::filter(ds_meta, study_label == "Schmiedel_2018", quant_method == "ge")

gtex_majiq_dsids = makeDatasetMetadata("Schmiedel_2018", gtex_tissues$sample_group, quant_methods = c("majiq"), 26, 1200)
write.table(gtex_majiq_dsids, "~/Downloads/new_dataset_id_map.tsv", sep = "\t", row.names = F, quote = F)

tissue_meta = dplyr::select(gtex_tissues, sample_group, tissue_id, tissue_label, condition_label, sample_size, pmid, study_type)

gtex_majiq_meta = dplyr::left_join(gtex_majiq_dsids, tissue_meta, by = "sample_group") %>%
  dplyr::select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size, quant_method, pmid, study_type)
write.table(gtex_majiq_meta, "~/Downloads/new_meta.tsv", sep = "\t", row.names = F, quote = F)



#Add MAJIQ ids to all other RNA-seq datasets
ds_meta_r7 = readr::read_tsv("~/projects/eQTL-Catalogue-resources/data_tables/dataset_metadata_r7.tsv")
ds_meta_upcoming = readr::read_tsv("~/projects/eQTL-Catalogue-resources/data_tables/dataset_metadata_upcoming.tsv")

#Quant methods with ge
ge_datasets = dplyr::filter(ds_meta_r7, quant_method == "ge", study_type == "bulk")

#Existing majiq studies
majiq_studies = ds_meta_upcoming %>% dplyr::select(study_id, quant_method) %>% dplyr::filter(quant_method == "majiq") %>% dplyr::distinct()

#missing majiq
no_majq_datasets = dplyr::anti_join(ge_datasets, majiq_studies, by = "study_id")

#Make new dataset ids
dataset_ids = makeDatasetMetadata("AAAA", no_majq_datasets$sample_group, quant_methods = "majiq", 
                                  study_index_start = 1, dataset_index_start = 1220)

new_majiq_datasets = dplyr::mutate(no_majq_datasets, dataset_id = dataset_ids$dataset_id, quant_method = "majiq")
write.table(new_majiq_datasets, "~/Downloads/new_meta.tsv", sep = "\t", row.names = F, quote = F)
