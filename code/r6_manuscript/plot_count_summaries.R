library("dplyr")

rnaseq_study_names <- c("Alasoo_2018", "BLUEPRINT_SE", "BLUEPRINT_PE", "BrainSeq", "GTEx", "FUSION", "GENCORD", "GEUVADIS", 
                        "HipSci", "Lepik_2017", "Nedelec_2016", "Quach_2016", "ROSMAP", "Schmiedel_2018", "Schwartzentruber_2018", 
                        "TwinsUK", "van_de_Bunt_2015", "Bossini-Castillo_2019", "CAP", "CommonMind", "Peng_2018", "PhLiPS", 
                        "iPSCORE","Braineac2","Steinberg_2020","Young_2019")
microarray_study_names <- c("CEDAR", "Fairfax_2014", "Kasela_2017", "Naranbhai_2015", "Fairfax_2012", "Gilchrist_2021")
file_names = c(rnaseq_study_names, microarray_study_names)
meta_path = "../SampleArcheology/studies/cleaned/"

old_studies = c("ROSMAP","BrainSeq","TwinsUK","FUSION","BLUEPRINT","Quach_2016",
                "Schmiedel_2018","GENCORD","GEUVADIS","Alasoo_2018",
                "Nedelec_2016","Lepik_2017","HipSci","van_de_Bunt_2015",
                "Schwartzentruber_2018","GTEx","CEDAR","Fairfax_2012",
                "Fairfax_2014","Kasela_2017","Naranbhai_2015")



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
  dplyr::mutate(dataset_name = paste0(study_file, "_", qtl_group)) %>% 
  dplyr::mutate(old_study = study %in% old_studies)

new_studies = setdiff(final_meta$study %>% unique(), old_studies)

final_meta_sum <- final_meta %>% 
  dplyr::group_by(old_study) %>% 
  dplyr::summarise(sample_count = length(unique(sample_id)),
                   donor_count = length(unique(genotype_id)),
                   dataset_count = length(unique(dataset_name)),
                   study_count = length(unique(study))) 

final_meta_sum <- final_meta_sum %>% 
  dplyr::bind_rows(final_meta_sum %>% colSums()) %>% 
  dplyr::mutate(state = c("new", "R3", "R6")) %>% 
  dplyr::filter(state != "new") %>% 
  dplyr::mutate(state = factor(state, levels = c("R3", "R6")))


library(ggplot2)
library(patchwork)

study_counts_plot = ggplot(final_meta_sum, aes(x=state, y=study_count)) +
  geom_bar(stat="identity", fill="#56B4E9")+
  geom_text(aes(label=study_count), vjust=1.6, color="white", size=3.5) +
  labs(x = "Study counts", y=NULL)+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

dataset_counts_plot = ggplot(final_meta_sum, aes(x=state, y=dataset_count)) +
  geom_bar(stat="identity", fill="#56B4E9")+
  geom_text(aes(label=dataset_count), vjust=1.6, color="white", size=3.5) +
  labs( x = "Dataset counts", y=NULL)+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

donor_counts_plot = ggplot(final_meta_sum, aes(x=state, y=donor_count)) +
  geom_bar(stat="identity", fill="#56B4E9")+
  geom_text(aes(label=donor_count), vjust=1.6, color="white", size=3.5) +
  labs( x = "Donor counts", y=NULL)+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

sample_counts_plot = ggplot(final_meta_sum, aes(x=state, y=sample_count)) +
  geom_bar(stat="identity", fill="#56B4E9")+
  geom_text(aes(label=sample_count), vjust=1.6, color="white", size=3.5) +
  labs(x = "Sample counts", y=NULL)+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

all_counts_plot = study_counts_plot | dataset_counts_plot | donor_counts_plot | sample_counts_plot

ggsave(filename = "manuscript_update/plots/counts_plot_nogrid_freeze_251122.pdf", plot = all_counts_plot, device = "pdf", width = 8, height = 4)
ggsave(filename = "manuscript_update/plots/counts_plot_nogrid_freeze_251122.png", plot = all_counts_plot, device = "png", width = 8, height = 4)


message("Number of studies: ", length(unique(final_meta$study)))
message("Number of datasets: ", length(unique(final_meta$dataset_name)))
message("Number of donors:", length(unique(final_meta$genotype_id)))

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
all_datasets = dplyr::bind_rows(rnaseq_datasets, microarray_datasets)

#Assing study and dataset ids
dataset_id_map = readr::read_tsv("data_tables/dataset_id_map.tsv")
dataset_metadata = dplyr::rename(all_datasets, study_label = study, 
                                 sample_group = qtl_group,
                                 tissue_id = tissue_ontology_id) %>% 
  dplyr::left_join(dataset_id_map, by = c("study_label", "quant_method", "sample_group")) %>%
  dplyr::select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size, quant_method) %>%
  dplyr::arrange(study_id, dataset_id)
write.table(dataset_metadata, "data_tables/dataset_metadata.tsv", sep = "\t", row.names = F, quote = F)


