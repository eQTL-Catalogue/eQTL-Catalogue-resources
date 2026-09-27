library("readr")
library("ggplot2")
library("dplyr")
library("data.table")

import_qtl_perm_res <- function(file_path, p_fdr_thresh = 0.01){
  col_types = "cciiccddddd"
  
  table <- readr::read_tsv(file_path, col_types = col_types) %>%
    dplyr::filter(!is.na(p_beta)) %>% 
    dplyr::mutate(txrev_event_type = gsub(pattern = ".*\\.", replacement = "", x = molecular_trait_object_id)) %>% 
    dplyr::group_by(txrev_event_type) %>% 
    dplyr::mutate(p_bonferroni = pvalue*n_variants*n_traits) %>%
    dplyr::mutate(p_bonferroni = pmin(p_bonferroni,1)) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr")) %>%
    dplyr::filter(p_fdr <= p_fdr_thresh) %>%
    dplyr::arrange(p_fdr) %>% 
    dplyr::ungroup()
  
  return(table)
}

V5_txrev_files_path <- "/Users/kerimov/Work/temp_files/V5_qtlmap_permutations"
V5_txrev_files = list.files(V5_txrev_files_path, pattern = "*txrev*", full.names = T)

V5_txrev_df = data.frame()
for (V5_txrev_file in V5_txrev_files) {
  message(" ## Reading file: ", V5_txrev_file)
  V5_txrev_sign <- import_qtl_perm_res(file_path = V5_txrev_file, p_fdr_thresh = 0.01)
  V5_txrev_sign <- V5_txrev_sign %>% 
    dplyr::mutate(txrev_event_type = gsub(pattern = ".*\\.", replacement = "", x = molecular_trait_object_id))
  
  study_group <- sub(pattern = V5_txrev_files_path, replacement = "", x = V5_txrev_file)
  study_group <- sub(pattern = ".permuted.tsv.gz", replacement = "", x = study_group)
  study_group <- sub(pattern = "_txrev_", replacement = ".", x = study_group)
  study_group <- sub(pattern = "/", replacement = "", x = study_group)
  
  V5_txrev_sign <- V5_txrev_sign %>% 
    dplyr::mutate(study_group = study_group)
  V5_txrev_df <- V5_txrev_df %>% rbind(V5_txrev_sign)
  
}

nrow(V5_txrev_df)
V5_txrev_df %>% colnames()
V5_txrev_df$txrev_event_type %>% table()
V5_txrev_df <- V5_txrev_df %>% 
  mutate(release = "V5")
#############################################################################
# 
# V6_txrev_files_path <- "/Users/kerimov/Work/temp_files/V6_qtlmap_permutations/txrev/"
# V6_txrev_files = list.files(V6_txrev_files_path, pattern = "*txrev*", full.names = T)
# 
# V6_txrev_df = data.frame()
# for (V6_txrev_file in V6_txrev_files) {
#   message(" ## Reading file: ", V6_txrev_file)
#   V6_txrev_sign <- import_qtl_perm_res(file_path = V6_txrev_file, p_fdr_thresh = 0.01)
#   V6_txrev_sign <- V6_txrev_sign %>% 
#     dplyr::mutate(txrev_event_type = gsub(pattern = ".*\\.", replacement = "", x = molecular_trait_object_id))
#   
#   study_group <- sub(pattern = V6_txrev_files_path, replacement = "", x = V6_txrev_file)
#   study_group <- sub(pattern = ".permuted.tsv.gz", replacement = "", x = study_group)
#   study_group <- sub(pattern = "_txrev_", replacement = ".", x = study_group)
#   study_group <- sub(pattern = "/", replacement = "", x = study_group)
#   
#   V6_txrev_sign <- V6_txrev_sign %>% 
#     dplyr::mutate(study_group = study_group)
#   V6_txrev_df <- V6_txrev_df %>% rbind(V6_txrev_sign)
#   
# }
# 
# nrow(V6_txrev_df)
# V6_txrev_df %>% colnames()
# V6_txrev_df$txrev_event_type %>% table()

master_data_txrev = master_data %>% 
  dplyr::mutate(study_group = paste0(study_label,".",sample_group)) %>% 
  dplyr::select(dataset_id, quant_method, study_group)

V6_txrev_df = V6_all_clean_df %>% 
  dplyr::left_join( master_data_txrev) %>% 
  dplyr::filter(quant_method == "txrev") %>% 
  dplyr::mutate(txrev_event_type = gsub(pattern = ".*\\.", replacement = "", x = molecular_trait_object_id)) %>% 
  dplyr::mutate(release = "V6")


V5_txrev_gene_df = V5_txrev_df %>% 
  mutate(gene_id = gsub(pattern = "\\..*", replacement = "", x = molecular_trait_object_id))

avail_colnames_both = intersect(colnames(V5_txrev_gene_df), colnames(V6_txrev_df))
txrev_df_all <- rbind(V5_txrev_gene_df %>% dplyr::select(all_of(avail_colnames_both)), 
                      V6_txrev_df %>% dplyr::select(all_of(avail_colnames_both)))

txrev_event_types = data.frame(txrev_event_type = c("contained", "downstream", "upstream"),
                               txrev_event_type_public = c("splicing", "3' end usage", "promoter usage"))

txrev_df_all_summ <- txrev_df_all %>% 
  dplyr::group_by(study_group, release, txrev_event_type) %>% 
  dplyr::summarise(.groups = "keep", e_genes_count = length(unique(gene_id))) %>% 
  dplyr::left_join(txrev_event_types)

txrev_df_all_summ %>% 
  dplyr::filter(txrev_event_type == "upstream") %>% 
  dplyr::group_by(release) %>% 
  dplyr::summarise(mean_egenes = mean(e_genes_count), sd_egenes = sd(e_genes_count)) 

txrev_event_type_egene_counts_plot = ggplot(txrev_df_all_summ, aes(x=txrev_event_type_public, y=e_genes_count, fill=release)) +
  geom_boxplot() + 
  labs( x = "Txrevise event type", y= "Count of unique eGenes per dataset") + #title = "Count of unique eGenes per dataset",
  theme_bw()

ggsave(txrev_event_type_egene_counts_plot, filename = "manuscript_update/plots/txrev_event_type_egene_counts_plot_freeze.pdf", device = "pdf", width = 7, height = 5)
ggsave(txrev_event_type_egene_counts_plot, filename = "manuscript_update/plots/txrev_event_type_egene_counts_plot_freeze.png", device = "png", width = 7, height = 5)
