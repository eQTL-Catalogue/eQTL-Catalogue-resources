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

p_fdr_thresh = 0.01
#############################################################################

master_data = readr::read_tsv("data_tables/dataset_metadata.tsv")

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


V6_ge_sign_counts = V6_all_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::filter(quant_method == "ge") %>% 
  dplyr::mutate(gene_id = molecular_trait_object_id)

V6_ma_sign_counts = V6_all_df %>% 
  dplyr::left_join(master_data %>% select(dataset_id, quant_method, sample_size, study_label)) %>% 
  dplyr::filter(quant_method == "microarray") %>% 
  dplyr::mutate(gene_id = molecular_trait_object_id)

no_x_chrom_studies = c("Nedelec_2016", "Schmiedel_2018","Peng_2018","van_de_Bunt_2015",
                       "Quach_2016","BrainSeq","CAP","TwinsUK","BLUEPRINT",
                       "GENCORD","Lepik_2017","ROSMAP","Steinberg_2020")

V6_ge_X_sign_df = V6_ge_sign_counts %>% 
  dplyr::bind_rows(V6_ma_sign_counts) %>% 
  dplyr::filter(chromosome == "X")

no_x_master_data = master_data %>% 
  dplyr::filter(quant_method %in% c("ge", "microarray")) %>% 
  dplyr::filter(study_label %in% no_x_chrom_studies) %>% 
  dplyr::select(study_label, dataset_id, sample_size) %>% 
  dplyr::mutate(e_genes_count = 0, has_x_chrom = FALSE )

ge_df_all_summ <- V6_ge_X_sign_df %>% 
  dplyr::group_by(study_label, dataset_id, sample_size) %>% 
  dplyr::summarise(.groups = "keep", e_genes_count = n()) %>% 
  dplyr::mutate(has_x_chrom = !study_label %in% no_x_chrom_studies) %>% 
  dplyr::bind_rows(no_x_master_data)

readr::write_tsv(ge_df_all_summ, file = "manuscript_update/X_chrom/x_chrom_egene_summ.tsv")

X_chrom_egene_count_plot <- ggplot(ge_df_all_summ, aes(x=sample_size, y=e_genes_count, color = has_x_chrom)) +
  geom_point(alpha = 0.7) +
  labs(x = "Sample size", y="No of unique genes in X cromosome", color = "Dataset contains X chromosome") +
  theme_bw()+
  scale_colour_manual(values=c("darkred", "darkblue")) +
  theme(legend.position = "none") +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave(filename = "manuscript_update/plots/Xchrom_egenes_nogrid_freeze_251122.pdf", plot = X_chrom_egene_count_plot, device = "pdf", height = 4, width = 4)
ggsave(filename = "manuscript_update/plots/Xchrom_egenes_nogrid_freeze_251122.png", plot = X_chrom_egene_count_plot, device = "png", height = 4, width = 4)

counts_and_x = (study_counts_plot | dataset_counts_plot | donor_counts_plot | sample_counts_plot | plot_spacer() | X_chrom_egene_count_plot) + plot_layout(nrow = 1, widths = c(1,1,1,1,0.1,2))  

# ggsave(filename = "manuscript_update/plots/x_chrom_egenes_count.pdf", plot = X_chrom_egene_count_plot, device = "pdf", height = 4, width = 4)
ggsave(filename = "manuscript_update/plots/counts_and_Xchrom_nogrid.pdf", plot = counts_and_x, device = "pdf", height = 4, width = 12)
ggsave(filename = "manuscript_update/plots/counts_and_Xchrom_nogrid.png", plot = counts_and_x, device = "png", height = 4, width = 12)
  
figure1 = (counts_and_x / counts_all_quants)+ plot_layout(nrow = 2, heights = c(3,2)) 
ggsave(filename = "manuscript_update/plots/Figure1.pdf", plot = figure1, device = "pdf", height = 7, width = 15)
ggsave(filename = "manuscript_update/plots/Figure1.png", plot = figure1, device = "png", height = 7, width = 15)

figure1_facet = (counts_and_x / all_egene_count_facet_plot) + plot_layout(nrow = 2, heights = c(3,2))  
ggsave(filename = "manuscript_update/plots/figure1_facet_nogrid.pdf", plot = figure1_facet, device = "pdf", height = 7, width = 13)
ggsave(filename = "manuscript_update/plots/figure1_facet_nogrid.png", plot = figure1_facet, device = "png", height = 7, width = 13)
