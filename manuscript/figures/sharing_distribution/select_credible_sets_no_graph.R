
library("dplyr")
library("readr")
library("ggplot2")

#set the required column names
colnames_susie <- c("phenotype_id", "variant_id","chr","pos","ref","alt", "cs_id",
                    "cs_index", "finemapped_region", "pip", "z", "cs_min_r2",
                    "cs_avg_r2", "cs_size", "posterior_mean", "posterior_sd", "cs_log10bf")

# read phenotype metadata and extract the needed columns
gene_pheno_meta <- read_tsv("../../temp_files/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz")
gene_pheno_meta <- gene_pheno_meta %>% dplyr::select(phenotype_id, gene_id)

# get the names of the needed susie files
gene_susie_files <- list.files(path = "../../temp_files/susie", pattern = "*ge.purity_filtered.txt.gz", full.names = T)

df_merged <- data.frame()

# for each dataset susie file
for (gene_susie_file in gene_susie_files) {
  # extract the name of the dataset from the file name
  dataset_name <- gsub(x = basename(gene_susie_file), pattern = "_ge.purity_filtered.txt.gz", replacement = "")
  
  # read the susie file
  susie_gene_df <- read_tsv(gene_susie_file, col_names = colnames_susie, skip = 1)
  
  # group the credible sets for the selected DATASET by gene_id and 
  # select the credible set of that gene with smallest credible set size and highest pip
  # Since we are selected from one dataset this selection is not biased for sample size
  susie_gene_df_selected_pip <- susie_gene_df %>% 
    left_join(gene_pheno_meta) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::arrange(cs_size, desc(pip)) %>% 
    dplyr::slice(1) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(dataset = dataset_name)
  
  # select the whole credible set (cs_size number of QTLs) which we selected 
  # per gene in the previous step
  susie_gene_cs_df <- susie_gene_df %>% 
    dplyr::filter(cs_id %in% susie_gene_df_selected_pip$cs_id) %>% 
    left_join(gene_pheno_meta) %>% 
    dplyr::mutate(dataset = dataset_name)
  
  # merge the selected credible sets of all of the datasets
  df_merged <- rbind(df_merged, susie_gene_cs_df)
}
nrow(df_merged)
ncol(df_merged)

# write_tsv(df_merged, "../../temp_files/eQTL_sharing/selected_crediblesets_gene_all_datasets_V8.tsv")

# generate cc_id column to have unique id for each dataset-credible_set
df_merged_cc_id <- df_merged %>% mutate(cc_id = paste0(dataset,"+",cs_id))

# At this point df_merged_cc_id can have different credible sets for the same gene but for different datasets
# we need only one credible set per gene 
# Since we can not choose a credible set for a gene across datasets (because sample size of datasets will bias the selection)
# we select one credible set per gene randomly across datasets
df_merged_random_select <- df_merged_cc_id %>%
  dplyr::group_by(gene_id) %>%
  dplyr::sample_n(1) %>%
  dplyr::ungroup()
#  23,248 genes and rows

df_merged_csuid_selected <- df_merged_cc_id %>% 
  filter(cc_id %in% df_merged_random_select$cc_id)
# 917,057 qtls from 23,248 credible sets

write_tsv(df_merged_csuid_selected, "../../temp_files/eQTL_sharing/ready_cc_id_V8/gene_cc_id_V8.tsv")
