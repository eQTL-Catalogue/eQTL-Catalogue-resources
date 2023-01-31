library("readr")
library("ggplot2")
library("dplyr")
library("data.table")

import_qtl_perm_res <- function(file_path){
  col_types = "cciiccddddd"
  table = readr::read_tsv(file_path, col_types = col_types) %>%
    dplyr::filter(!is.na(p_beta)) %>%
    dplyr::mutate(p_bonferroni = pvalue*n_variants*n_traits) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::arrange(p_fdr)
  return(table)
}

old_studies = c("ROSMAP","BrainSeq","TwinsUK","FUSION","BLUEPRINT","Quach_2016",
                "Schmiedel_2018","GENCORD","GEUVADIS","Alasoo_2018",
                "Nedelec_2016","Lepik_2017","HipSci","van_de_Bunt_2015",
                "Schwartzentruber_2018","GTEx","CEDAR","Fairfax_2012",
                "Fairfax_2014","Kasela_2017","Naranbhai_2015")



p_fdr_thresh = 0.01
#############################################################################
# sample_sizes = readr::read_tsv("/Users/kerimov/Work/GitHub/eQTL-Catalogue-resources/data_tables/V6/sample_size_by_dataset.tsv") %>% 
#   dplyr::mutate(study_group = paste0(study, ".", qtl_group)) %>% 
#   dplyr::select(study_group, study, sample_size) %>% 
#   dplyr::mutate(study = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study)) %>%
#   dplyr::mutate(study = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study)) %>%
#   dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study_group)) %>% 
#   dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study_group))

master_data = readr::read_tsv("data_tables/dataset_metadata.tsv")


# quant_methods = c("ge", "exon", "tx", "txrev", "leafcutter")
V6_perm_files_path <- "/Users/kerimov/Work/temp_files/__freeze_251122_/"
V6_perm_files = list.files(V6_perm_files_path, full.names = T)

V6_all_df = data.frame()
for (V6_perm_file in V6_perm_files) {
  message(" ## Reading file: ", V6_perm_file)
  V6_sign <- import_qtl_perm_res(file_path = V6_perm_file)
  V6_sign <- V6_sign %>%
    dplyr::filter(p_fdr <= p_fdr_thresh)

  study_group <- sub(pattern = V6_perm_files_path, replacement = "", x = V6_perm_file)
  study_group <- sub(pattern = ".permuted.tsv.gz", replacement = "", x = study_group)
  # study_group <- sub(pattern = paste0("_",quant_method,"_"), replacement = ".", x = study_group)
  study_group <- sub(pattern = "/", replacement = "", x = study_group)

  V6_sign <- V6_sign %>%
    dplyr::mutate(dataset_id = study_group) 

  V6_all_df <- V6_all_df %>% 
    dplyr::bind_rows(V6_sign)
}

# 
# V6_all_df = V6_all_df %>% 
#   dplyr::mutate(study = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study)) %>% 
#   dplyr::mutate(study = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study)) %>% 
#   dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study_group)) %>% 
#   dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study_group)) 

V6_all_clean_df <- V6_all_df %>%
  dplyr::mutate(gene_id = gsub(pattern = ".contained|.upstream|.downstream", replacement = "", x = molecular_trait_object_id))

# Plot GE e_genes
V6_ge_df_egene_counts = V6_all_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::filter(quant_method %in% c("ge", "microarray")) %>% 
  dplyr::group_by(quant_method, study_label, dataset_id, sample_size) %>% 
  dplyr::summarise(unique_egene_count = length(unique(molecular_trait_id))) %>% 
  # dplyr::left_join(sample_sizes) %>% 
  dplyr::mutate(old_study = study_label %in% old_studies) %>% 
  dplyr::mutate(quant_method = "ge")

ge_egene_count_plot <- ggplot(V6_ge_df_egene_counts, aes(x=sample_size, y=unique_egene_count, color = old_study)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="Gene expression\neGene counts", color = "State") +
  theme_bw()+
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  # scale_f(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none") + 
  ylim(0, 11000) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())



# Plot exon e_genes
V6_exon_df_egene_counts = V6_all_clean_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::filter(quant_method == "exon") %>% 
  dplyr::group_by(quant_method, study_label, dataset_id, sample_size) %>% 
  dplyr::summarise(unique_egene_count = length(unique(molecular_trait_id))) %>% 
  dplyr::mutate(old_study = study_label %in% old_studies)

exon_egene_count_plot <- ggplot(V6_exon_df_egene_counts, aes(x=sample_size, y=unique_egene_count, color = old_study)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="Exon expression\neGene counts", color = "State") +
  theme_bw()+
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none") +
  ylim(0, 11000)+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# Plot TX e_genes
V6_tx_df_egene_counts = V6_all_clean_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::filter(quant_method == "tx") %>% 
  dplyr::group_by(quant_method, study_label, dataset_id, sample_size) %>% 
  dplyr::summarise(unique_egene_count = length(unique(molecular_trait_id))) %>% 
  dplyr::mutate(old_study = study_label %in% old_studies)

tx_egene_count_plot <- ggplot(V6_tx_df_egene_counts, aes(x=sample_size, y=unique_egene_count, color = old_study)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="Transcript usage\neGene counts", color = "State") +
  theme_bw()+
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none") + 
  ylim(0, 11000)+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# Plot TXREV e_genes
V6_txrev_df_egene_counts = V6_all_clean_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::filter(quant_method == "txrev") %>% 
  dplyr::group_by(quant_method, study_label, dataset_id, sample_size) %>% 
  dplyr::summarise(unique_egene_count = length(unique(molecular_trait_id))) %>% 
  dplyr::mutate(old_study = study_label %in% old_studies)

txrev_egene_count_plot <- ggplot(V6_txrev_df_egene_counts, aes(x=sample_size, y=unique_egene_count, color = old_study)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="TxRev usage\neGene counts", color = "State") +
  theme_bw()+
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none")+ 
  ylim(0, 11000)+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


# V6_leafcutter_perm_files = list.files("/Users/kerimov/Work/temp_files/V6_qtlmap_permutations/leafcutter", full.names = T)
leafcutter_meta_map = readr::read_tsv("/Users/kerimov/Work/temp_files/V6_qtlmap_permutations/leafcutter_meta_map.tsv") %>% 
  dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study_group)) %>%
  dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study_group))

V6_lc_df = data.frame()
# prepare leafcutter eGenes

lc_master_data = master_data %>% 
  dplyr::filter(quant_method == "leafcutter") %>% 
  dplyr::mutate(qtl_group = paste0(study_label, ".", sample_group))

for (index in 1:nrow(lc_master_data)) {
  message(" ## Reading file: ", lc_master_data$qtl_group[index])
  V6_sign <- V6_all_clean_df %>% 
    dplyr::filter(dataset_id == lc_master_data$dataset_id[index])
  
  leafcutter_meta_file = leafcutter_meta_map %>% 
    dplyr::filter(study_group == lc_master_data$qtl_group[index]) %>% 
    dplyr::pull(metadata_file)
  
  leafcutter_meta = readr::read_tsv(leafcutter_meta_file, col_types = "cccccddiccidddddddd") %>% 
    dplyr::select(molecular_trait_id = phenotype_id, gene_id)
  
  V6_sign = V6_sign %>% 
    dplyr::select(-gene_id) %>% 
    dplyr::left_join(leafcutter_meta, by = "molecular_trait_id")
  
  V6_lc_df <- V6_lc_df %>% 
    dplyr::bind_rows(V6_sign)
}

# V6_quant_method_df = V6_quant_method_df %>%
  # dplyr::mutate(study = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study)) %>%
  # dplyr::mutate(study = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study)) %>%
  # dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT", x = study_group)) %>%
  # dplyr::mutate(study_group = gsub(pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT", x = study_group))

# Plot LEAFCUTTER e_genes
V6_leafcutter_df_egene_counts = V6_lc_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::group_by(quant_method, study_label, dataset_id, sample_size) %>% 
  dplyr::summarise(unique_egene_count = length(unique(molecular_trait_id))) %>% 
  dplyr::mutate(old_study = study_label %in% old_studies) 

leafcuuter_egene_count_plot <- ggplot(V6_leafcutter_df_egene_counts, aes(x=sample_size, y=unique_egene_count, color = old_study)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="Splicing usage (Leafcutter)\neGene counts", color = "State") +
  theme_bw()+
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none") + 
  ylim(0, 11000)+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

all_egene_counts = V6_ge_df_egene_counts %>% 
  dplyr::bind_rows(V6_exon_df_egene_counts) %>% 
  dplyr::bind_rows(V6_tx_df_egene_counts) %>% 
  dplyr::bind_rows(V6_txrev_df_egene_counts) %>% 
  dplyr::bind_rows(V6_leafcutter_df_egene_counts)

all_egene_count_facet_plot <- ggplot(all_egene_counts, aes(x=sample_size, y=unique_egene_count, color = old_study)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="No of unique genes", color = "State") +
  theme_bw()+
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  theme(legend.position = "none") + 
  ylim(0, 11000) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  facet_grid(cols = vars(factor(quant_method, levels = c("ge", "exon", "tx", "txrev", "leafcutter")))) +
  theme(strip.background =element_rect(fill="white", color = "white"))

ggsave(filename = "manuscript_update/plots/all_egene_count_facet_plot_nogrid_freeze_251122.pdf", plot = all_egene_count_facet_plot, device = "pdf", height = 3, width = 12)
ggsave(filename = "manuscript_update/plots/all_egene_count_facet_plot_nogrid_freeze_251122.png", plot = all_egene_count_facet_plot, device = "png", height = 3, width = 12)


library(patchwork)
counts_all_quants= ge_egene_count_plot | exon_egene_count_plot | leafcuuter_egene_count_plot | tx_egene_count_plot | txrev_egene_count_plot # + plot_layout(nrow = 1, widths = c(1,1,1,1,0.2,2))

# ggsave(filename = "manuscript_update/plots/x_chrom_egenes_count.pdf", plot = X_chrom_egene_count_plot, device = "pdf", height = 4, width = 4)
ggsave(filename = "manuscript_update/plots/counts_all_quants_ylim_freeze_251122.pdf", plot = counts_all_quants, device = "pdf", height = 4, width = 20)
ggsave(filename = "manuscript_update/plots/counts_all_quants_ylim_freeze_251122.png", plot = counts_all_quants, device = "png", height = 4, width = 20)
