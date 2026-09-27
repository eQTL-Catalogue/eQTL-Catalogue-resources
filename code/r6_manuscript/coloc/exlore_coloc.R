library("readr")
library("ggplot2")
library("dplyr")
library("data.table")

GWAS_traits = c("UKBB.VitD","UKBB.Mono","UKBB.Lym","UKBB.Height","UKBB.WBC","UKBB.TG","UKBB.TC","UKBB.RBC","UKBB.Plt","UKBB.MCV","UKBB.MCP","UKBB.LipoA","UKBB.LDLC","UKBB.HDLC","UKBB.Eosino","UKBB.CAD","UKBB.BMI")

coloc_path = "/Users/kerimov/Work/temp_files/coloc_v5_results_merged/"
susie_path = "/Users/kerimov/Work/temp_files/susie_all"
datasets = list.files(path = coloc_path)
# FUSION_leafcutter_adipose_naive_UKBB.LDLC.coloc.v5 = readr::read_tsv("/Users/kerimov/Work/temp_files/coloc_results_Ralf/FUSION_leafcutter_adipose_naive_UKBB.LDLC.coloc.v5.txt.gz") %>% 
#   dplyr::filter(PP.H4.abf > 0.8)
# FUSION_leafcutter_adipose_naive_UKBB.LDLC.coloc.v5$hit1 

all_coverage_plots = readr::read_tsv("/Users/kerimov/Work/temp_files/coloc_results_Ralf/all_leacutter_coverage_plots.tsv", col_names = "plot_path")

all_coverage_plots_mut = all_coverage_plots %>% 
  dplyr::mutate(plot_path = gsub(pattern = "./", replacement = "", x = plot_path, fixed = T)) %>% 
  dplyr::mutate(dataset = str_split_fixed(string = plot_path, pattern = "/", n = 2)[,1]) %>% 
  dplyr::filter(!dataset %in% c("logs","pipeline_info")) %>% 
  dplyr::mutate(signal_name = str_split_fixed(string = plot_path, pattern = "/", n = 2)[,2]) %>% 
  dplyr::mutate(molecular_trait_id = str_split_fixed(string = signal_name, pattern = "&", n = 3)[,1]) %>% 
  dplyr::mutate(variant = str_split_fixed(string = signal_name, pattern = "&", n = 3)[,2]) %>% 
  dplyr::mutate(gene = str_split_fixed(string = signal_name, pattern = "&", n = 3)[,3]) %>% 
  dplyr::mutate(dataset_only = gsub(pattern = "_leafcutter", replacement = "", x = dataset))

coverage_plots_exist_only = data.frame()

# , 
# col_names = c("gene1","gene2","hit1","hit2","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","idx1","idx2"),
# col_types = "ccccdddddddd"

for (dataset in datasets) {
  message(" ## Processing: ", dataset)
  coloc_file_path = file.path(coloc_path, dataset)
  coloc_files_for_dataset = list.files(coloc_file_path)
  sign_coloc_signals = data.frame()
  for (coloc_file_for_dataset in coloc_files_for_dataset) {
    sign_coloc_signals_temp = readr::read_tsv(file.path(coloc_file_path, coloc_file_for_dataset)) 
    
    if (nrow(sign_coloc_signals_temp) == 0) {
      next
    }
    sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
      dplyr::filter(PP.H4.abf > 0.8)
    
    gwas_trait = gsub(pattern = dataset,replacement = "", coloc_file_for_dataset)
    gwas_trait = gsub(pattern = ".coloc.v5.txt.gz", replacement = "", gwas_trait)
    gwas_trait = gsub(pattern = "_", replacement = "", gwas_trait)
    sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
      dplyr::mutate(GWAS_trait = gwas_trait) %>% 
      dplyr::mutate(dataset = dataset)
    
    sign_coloc_signals = sign_coloc_signals %>% 
      dplyr::bind_rows(sign_coloc_signals_temp)
  }
  
  if (nrow(sign_coloc_signals) == 0) {
    next
  }
  susie_for_dataset = readr::read_tsv(paste0(file.path(susie_path, dataset), ".purity_filtered.txt.gz"), 
                                                          col_types = "cccdcccccdddddddd")
  
  susie_for_dataset_merged = susie_for_dataset %>% 
    dplyr::rename(gene1 = molecular_trait_id, hit1 = variant) %>% 
    dplyr::right_join(sign_coloc_signals, by = c("gene1", "hit1")) %>% 
    dplyr::rename(molecular_trait_id = gene1, variant = hit1) %>% 
    dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = molecular_trait_id))
  
  dataset_clean = gsub(pattern = "_leafcutter_", replacement = "_", x = dataset)
  
  all_coverage_plots_mut_dataset = all_coverage_plots_mut %>% 
    dplyr::filter(dataset_only == dataset_clean) 
  
  susie_for_dataset_with_cov_plot = susie_for_dataset_merged %>% 
    dplyr::select(-dataset) %>% 
    dplyr::left_join(all_coverage_plots_mut_dataset, by = c("molecular_trait_id", "variant"))
  
  message(" ## Writing : ", dataset_clean)
  readr::write_tsv(susie_for_dataset_with_cov_plot, file = paste0("colocs_to_covplots_output/", dataset_clean, ".tsv"))
  susie_for_dataset_with_cov_plot = susie_for_dataset_with_cov_plot %>% 
    dplyr::filter(!is.na(plot_path))
  
  readr::write_tsv(susie_for_dataset_with_cov_plot, file = paste0("colocs_to_covplots_output/", dataset_clean, "_exist_only.tsv"))
    
  coverage_plots_exist_only = coverage_plots_exist_only %>% 
    dplyr::bind_rows(susie_for_dataset_with_cov_plot)
  # intersect(all_coverage_plots_mut_dataset %>% colnames(), susie_for_dataset_merged %>% colnames())
}

readr::write_tsv(coverage_plots_exist_only, file = "coverage_plots_exist_only.tsv")

coverage_plots_exist_only_mut = coverage_plots_exist_only %>% 
  dplyr::mutate(full_path = paste0("/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_leafcutter_all/", plot_path)) %>% 
  dplyr::mutate(command = paste0("rsync -r '", full_path, "' ", "/gpfs/space/projects/eQTLCatalogue/coverage_plots/leafcutter_coloc_coverage_plots/", dataset, "/"))

write(coverage_plots_exist_only_mut$command, file = "commands_to_copy_covplots.sh")
# molecular_trait_id
# variant
# chromosome
# position
# ref
# alt
# cs_id
# cs_index
# finemapped_region
# pip
# z
# cs_min_r2
# cs_avg_r2
# cs_size
# posterior_mean
# posterior_sd
# cs_log10bf
FUSION_leafcutter_adipose_naive_susie = readr::read_tsv("/Users/kerimov/Work/temp_files/susie_all/FUSION_leafcutter_adipose_naive.purity_filtered.txt.gz", 
                                                        col_types = "cccdcccccdddddddd")

FUSION_leafcutter_adipose_naive_susie_merged = FUSION_leafcutter_adipose_naive_susie %>% 
  dplyr::rename(gene1 = molecular_trait_id, hit1 = variant) %>% 
  dplyr::right_join(FUSION_leafcutter_adipose_naive_UKBB.LDLC.coloc.v5) %>% 
  dplyr::rename(molecular_trait_id = gene1, variant = hit1) %>% 
  dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = molecular_trait_id))

available_coverage_plots = readr::read_tsv("/Users/kerimov/Work/temp_files/coloc_results_Ralf/FUSION_adipose_naive_leafcutter_plots.tsv", col_names = "signal_name")

available_coverage_plots_mut = available_coverage_plots %>% 
  dplyr::mutate(molecular_trait_id = str_split_fixed(string = signal_name, pattern = "&", n = 3)[,1]) %>% 
  dplyr::mutate(variant = str_split_fixed(string = signal_name, pattern = "&", n = 3)[,2]) %>% 
  dplyr::mutate(gene = str_split_fixed(string = signal_name, pattern = "&", n = 3)[,3])

available_coverage_plots_needed = FUSION_leafcutter_adipose_naive_susie_merged %>% 
  dplyr::left_join(available_coverage_plots_mut)

write_tsv(available_coverage_plots_needed, file = "available_coverage_plots_needed.tsv")


##############################################################################
## Read metadata files
metadata_files_df = readr::read_tsv("/Users/kerimov/Work/temp_files/leafcutter_metadata_109/metadata_dataset_mapping.tsv") %>% 
  dplyr::mutate(dataset  = paste0(study, "_", qtlset))

all_leafcutter_metadata = data.frame()
for (i in 1:nrow(metadata_files_df)) {
  message(" ## Processing: ", metadata_files_df$dataset[i])
  metadata_temp = readr::read_tsv(metadata_files_df$leafcutter_metadata[i], 
                                  col_types = "cccccddiccidddddddd")
  metadata_temp = metadata_temp %>% 
    dplyr::mutate(dataset = metadata_files_df$dataset[i]) %>% 
    dplyr::mutate(study = metadata_files_df$study[i])
  
  all_leafcutter_metadata = all_leafcutter_metadata %>%   
    dplyr::bind_rows(metadata_temp)
}


# gzfile = gzfile("data/leafcutter_metadata_all.tsv.gz", "w")
# write.table(x = all_leafcutter_metadata, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
# close(gzfile)

all_leafcutter_metadata$quant = "leafcutter"

all_leafcutter_metadata_needed = all_leafcutter_metadata %>% 
  dplyr::select(quant, dataset, gene1 = phenotype_id, gene_id, gene_name)

ge_pheno_meta = readr::read_tsv("/Users/kerimov/Work/temp_files/phenotype_metadata/gene_counts_Ensembl_105_phenotype_metadata.tsv.gz", 
                                col_types = "cccccddicciddi") %>% 
  dplyr::mutate(quant = "ge") %>% 
  dplyr::select(gene1 = phenotype_id, gene_id, gene_name)

exon_pheno_meta = readr::read_tsv("/Users/kerimov/Work/temp_files/phenotype_metadata/exon_counts_Ensembl_105_phenotype_metadata.tsv.gz", 
                                  col_types = "cccccddicciddi") %>% 
  dplyr::mutate(quant = "exon")%>% 
  dplyr::select(gene1 = phenotype_id, gene_id, gene_name)

tx_pheno_meta = readr::read_tsv("/Users/kerimov/Work/temp_files/phenotype_metadata/transcript_usage_Ensembl_105_phenotype_metadata.tsv.gz", 
                                col_types = "cccccddiccid") %>% 
  dplyr::mutate(quant = "tx") %>% 
  dplyr::select(gene1 = phenotype_id, gene_id, gene_name)

txrev_pheno_meta = readr::read_tsv("/Users/kerimov/Work/temp_files/phenotype_metadata/txrevise_Ensembl_105_phenotype_metadata.tsv.gz", 
                                   col_types = "cccccddiccid") %>% 
  dplyr::mutate(quant = "txrev")%>% 
  dplyr::select(gene1 = phenotype_id, gene_id, gene_name)
# -----------------------------------------------------------------------------

# Extract Vitamin D colocalisations 
GWAS_trait_oi = "UKBB.VitD"
leafcutter_colocs_path = "/Users/kerimov/Work/temp_files/coloc_results_all/results_leafcutter_UKBB_cc/coloc_v5_results_merged/"
leafcutter_coloc_files = list.files(path = leafcutter_colocs_path, pattern = paste0("*", GWAS_trait_oi, ".coloc.v5.txt.gz"), full.names = T, recursive = T)
leafcutter_colocs_df = data.frame(file_path = leafcutter_coloc_files) %>% 
  dplyr::mutate(filename = basename(file_path)) %>% 
  dplyr::mutate(dataset_quant = gsub(pattern = paste0("_", GWAS_trait_oi, ".coloc.v5.txt.gz"), replacement = "", x = filename)) %>% 
  dplyr::mutate(dataset = gsub(pattern = "_leafcutter_", replacement = "_", x = dataset_quant))

coloc_sign_threshold = 0.9

sign_coloc_signals = data.frame()
for (index in 1:nrow(leafcutter_colocs_df)) {
  message("# Reading: ", leafcutter_colocs_df[index, "dataset"], ", Count: ", index)
  sign_coloc_signals_temp = readr::read_tsv(leafcutter_colocs_df[index,"file_path"], col_types = "ccccdddddddd") 
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::filter(PP.H4.abf > coloc_sign_threshold) %>% 
    dplyr::mutate(dataset = leafcutter_colocs_df[index, "dataset"]) %>% 
    dplyr::mutate(quant = "leafcutter") 
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::left_join(all_leafcutter_metadata_needed, by = c("quant", "dataset", "gene1"))
  
  sign_coloc_signals <- sign_coloc_signals %>% 
    dplyr::bind_rows(sign_coloc_signals_temp)
}

#### Count unique GWAS loci 
count_unique_gene2_df = data.frame()
for (index in 1:nrow(leafcutter_colocs_df)) {
  message("# Reading: ", leafcutter_colocs_df[index, "dataset"])
  sign_coloc_signals_temp = readr::read_tsv(leafcutter_colocs_df[index,"file_path"], col_types = "ccccdddddddd")

  count_unique_gene2_df = count_unique_gene2_df %>%
    dplyr::bind_rows(sign_coloc_signals_temp %>% dplyr::mutate(dataset = leafcutter_colocs_df[index, "dataset"]))
}

tx_colocs_path = "/Users/kerimov/Work/temp_files/coloc_results_all/results_tx_UKBB_cc/coloc_v5_results_merged/"
tx_coloc_files = list.files(path = tx_colocs_path, pattern = paste0("*", GWAS_trait_oi, ".coloc.v5.txt.gz"), full.names = T, recursive = T)
tx_colocs_df = data.frame(file_path = tx_coloc_files) %>% 
  dplyr::mutate(filename = basename(file_path)) %>% 
  dplyr::mutate(dataset_quant = gsub(pattern = paste0("_", GWAS_trait_oi, ".coloc.v5.txt.gz"), replacement = "", x = filename)) %>% 
  dplyr::mutate(dataset = gsub(pattern = "_tx_", replacement = "_", x = dataset_quant))

for (index in 1:nrow(tx_colocs_df)) {
  message("# Reading: ", tx_colocs_df[index, "dataset"], ", Count: ", index)
  sign_coloc_signals_temp = readr::read_tsv(tx_colocs_df[index,"file_path"], col_types = "ccccdddddddd") 
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::filter(PP.H4.abf > coloc_sign_threshold) %>% 
    dplyr::mutate(dataset = tx_colocs_df[index, "dataset"]) %>% 
    dplyr::mutate(quant = "tx")
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::left_join(tx_pheno_meta, by = c("gene1"))
  
  sign_coloc_signals <- sign_coloc_signals %>% 
    dplyr::bind_rows(sign_coloc_signals_temp)
}

txrev_colocs_path = "/Users/kerimov/Work/temp_files/coloc_results_all/results_txrev_UKBB_cc/coloc_v5_results_merged/"
txrev_coloc_files = list.files(path = txrev_colocs_path, pattern = paste0("*", GWAS_trait_oi, ".coloc.v5.txt.gz"), full.names = T, recursive = T)
txrev_colocs_df = data.frame(file_path = txrev_coloc_files) %>% 
  dplyr::mutate(filename = basename(file_path)) %>% 
  dplyr::mutate(dataset_quant = gsub(pattern = paste0("_", GWAS_trait_oi, ".coloc.v5.txt.gz"), replacement = "", x = filename)) %>% 
  dplyr::mutate(dataset = gsub(pattern = "_txrev_", replacement = "_", x = dataset_quant))

for (index in 1:nrow(txrev_colocs_df)) {
  message("# Reading: ", txrev_colocs_df[index, "dataset"], ", Count: ", index)
  sign_coloc_signals_temp = readr::read_tsv(txrev_colocs_df[index,"file_path"], col_types = "ccccdddddddd") 
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::filter(PP.H4.abf > coloc_sign_threshold) %>% 
    dplyr::mutate(dataset = txrev_colocs_df[index, "dataset"]) %>% 
    dplyr::mutate(quant = "txrev")
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::left_join(txrev_pheno_meta, by = c("gene1"))
  
  sign_coloc_signals <- sign_coloc_signals %>% 
    dplyr::bind_rows(sign_coloc_signals_temp)
}


exon_colocs_path = "/Users/kerimov/Work/temp_files/coloc_results_all/results_exon_UKBB_cc/coloc_v5_results_merged/"
exon_coloc_files = list.files(path = exon_colocs_path, pattern = paste0("*", GWAS_trait_oi, ".coloc.v5.txt.gz"), full.names = T, recursive = T)
exon_colocs_df = data.frame(file_path = exon_coloc_files) %>% 
  dplyr::mutate(filename = basename(file_path)) %>% 
  dplyr::mutate(dataset_quant = gsub(pattern = paste0("_", GWAS_trait_oi, ".coloc.v5.txt.gz"), replacement = "", x = filename)) %>% 
  dplyr::mutate(dataset = gsub(pattern = "_exon_", replacement = "_", x = dataset_quant))

for (index in 1:nrow(exon_colocs_df)) {
  message("# Reading: ", exon_colocs_df[index, "dataset"], ", Count: ", index)
  sign_coloc_signals_temp = readr::read_tsv(exon_colocs_df[index,"file_path"], col_types = "ccccdddddddd") 
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::filter(PP.H4.abf > coloc_sign_threshold) %>% 
    dplyr::mutate(dataset = exon_colocs_df[index, "dataset"]) %>% 
    dplyr::mutate(quant = "exon")
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::left_join(exon_pheno_meta, by = c("gene1"))
  
  sign_coloc_signals <- sign_coloc_signals %>% 
    dplyr::bind_rows(sign_coloc_signals_temp)
}

ge_colocs_path = "/Users/kerimov/Work/temp_files/coloc_results_all/results_ge_UKBB_cc/coloc_v5_results_merged/"
ge_coloc_files = list.files(path = ge_colocs_path, pattern = paste0("*", GWAS_trait_oi, ".coloc.v5.txt.gz"), full.names = T, recursive = T)
ge_colocs_df = data.frame(file_path = ge_coloc_files) %>% 
  dplyr::mutate(filename = basename(file_path)) %>% 
  dplyr::mutate(dataset_quant = gsub(pattern = paste0("_", GWAS_trait_oi, ".coloc.v5.txt.gz"), replacement = "", x = filename)) %>% 
  dplyr::mutate(dataset = gsub(pattern = "_ge_", replacement = "_", x = dataset_quant)) %>% 
  dplyr::filter(stringr::str_detect(string = dataset, pattern = "microarray", negate = T))

for (index in 1:nrow(ge_colocs_df)) {
  message("# Reading: ", ge_colocs_df[index, "dataset"], ", Count: ", index)
  sign_coloc_signals_temp = readr::read_tsv(ge_colocs_df[index,"file_path"], col_types = "ccccdddddddd") 
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::filter(PP.H4.abf > coloc_sign_threshold) %>% 
    dplyr::mutate(dataset = ge_colocs_df[index, "dataset"]) %>% 
    dplyr::mutate(quant = "ge")
  
  if (nrow(sign_coloc_signals_temp) == 0) {
    next
  }
  
  sign_coloc_signals_temp = sign_coloc_signals_temp %>% 
    dplyr::left_join(ge_pheno_meta, by = c("gene1"))
  
  sign_coloc_signals <- sign_coloc_signals %>% 
    dplyr::bind_rows(sign_coloc_signals_temp)
}

readr::write_tsv(sign_coloc_signals, file = "data/VitD_sign_colocs_PP4_0.9.tsv")

sign_coloc_signals = readr::read_tsv("data/VitD_sign_colocs_PP4_0.9.tsv")

sign_coloc_signals_summary = sign_coloc_signals %>% 
  dplyr::group_by(gene2, idx2, gene_id, gene_name) %>%
  dplyr::mutate(leafcutter = "leafcutter" %in% unique(quant)) %>% 
  dplyr::mutate(ge = "ge" %in% unique(quant)) %>% 
  dplyr::mutate(tx = "tx" %in% unique(quant)) %>% 
  dplyr::mutate(txrev = "txrev" %in% unique(quant)) %>% 
  dplyr::mutate(exon = "exon" %in% unique(quant)) %>% 
  dplyr::summarise(ge, exon, tx, txrev, leafcutter, count = n()) %>% 
  dplyr::distinct()

readr::write_tsv(sign_coloc_signals_summary, file = "data/VitD_sign_colocs_PP4_0.9_summary_by_quant.tsv")

sign_coloc_signals_summary_3 = sign_coloc_signals %>% 
  dplyr::group_by(gene2, idx2) %>%
  dplyr::mutate(leafcutter = "leafcutter" %in% unique(quant)) %>% 
  dplyr::mutate(ge = "ge" %in% unique(quant)) %>% 
  dplyr::mutate(tx = "tx" %in% unique(quant)) %>% 
  dplyr::mutate(txrev = "txrev" %in% unique(quant)) %>% 
  dplyr::mutate(exon = "exon" %in% unique(quant)) %>% 
  dplyr::summarise(variant_hit2 = paste(unique(hit2), collapse = ","), gene_ids = paste(unique(gene_id), collapse = ","), gene_names = paste(unique(gene_name), collapse = ","), ge, exon, tx, txrev, leafcutter, count = n()) %>% 
  dplyr::distinct()

readr::write_tsv(sign_coloc_signals_summary_3, file = "data/VitD_sign_colocs_PP4_0.9_summary_by_quant_3.tsv")  

sign_coloc_signals_summary = sign_coloc_signals %>% 
  dplyr::group_by(gene2) %>% 
  dplyr::summarise(detected_by_quants = paste(unique(quant), collapse = ", "), dataset_count = length(unique(dataset)))
readr::write_tsv(sign_coloc_signals_summary, file = "data/VitD_sign_colocs_PP4_0.9_summary.tsv")

sign_coloc_signals_summary_by_dataset = sign_coloc_signals %>% 
  dplyr::group_by(gene2, dataset) %>% 
  dplyr::summarise(detected_by_quants = paste(unique(quant), collapse = ", "), dataset_count = length(unique(dataset)))
readr::write_tsv(sign_coloc_signals_summary_by_dataset, file = "data/VitD_sign_colocs_PP4_0.9_summary_by_dataset.tsv")


sign_coloc_signals_leafcutter = sign_coloc_signals %>% 
  dplyr::filter(quant == "leafcutter") %>% 
  dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = gene1, fixed = T)) %>% 
  # dplyr::mutate(coverage_plot_signal_start = paste0(molecular_trait_id, "&", hit1)) %>% 
  dplyr::mutate(dataset = paste0(dataset, "_", quant)) %>% 
  dplyr::left_join(all_coverage_plots_mut, by = c("dataset", "molecular_trait_id")) %>% 
  dplyr::mutate(full_path = paste0("/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_leafcutter_all/", plot_path)) %>% 
  dplyr::mutate(command = paste0("rsync -r '", full_path, "' ", "'/gpfs/space/projects/eQTLCatalogue/coverage_plots/VitD_leafcutter_coloc_coverage_plots_PP4_0.9/", gene2,"_idx", idx2, "/'"))

readr::write_tsv(sign_coloc_signals_leafcutter, file = "data/VitD_sign_colocs_PP4_0.9_coverage_plots.tsv")

write(sign_coloc_signals_leafcutter$command, file = "commands_to_copy_leafcutter_covplots.sh")

uniq_gene2_idx_combs = sign_coloc_signals %>% 
  dplyr::mutate(folder_name = paste0(gene2, "_idx", idx2)) %>% 
  dplyr::pull(folder_name) %>% 
  unique()

commands_make_gene2_folders = paste0("mkdir '/gpfs/space/projects/eQTLCatalogue/coverage_plots/VitD_leafcutter_coloc_coverage_plots_PP4_0.9/",uniq_gene2_idx_combs, "'")
write(commands_make_gene2_folders, file = "mkdir_VitD_GWAS_loci.sh")

sign_coloc_signals_exon = sign_coloc_signals %>% 
  dplyr::filter(quant == "exon") %>%
  dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = gene1, fixed = T)) %>% 
  # dplyr::mutate(coverage_plot_signal_start = paste0(molecular_trait_id, "&", hit1)) %>% 
  dplyr::mutate(dataset = paste0(dataset, "_", quant)) %>% 
  dplyr::mutate(copy_command = paste0("rsync -r '/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_exon_all/",dataset, "/",molecular_trait_id, "___", hit1,"' ", gene2,"_idx", idx2, "/", molecular_trait_id, "_", dplyr::row_number()))

write(sign_coloc_signals_exon$copy_command, file = "copy_exon_coverage_plots.sh")

sign_coloc_signals_tx = sign_coloc_signals %>% 
  dplyr::filter(quant == "tx") %>%
  dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = gene1, fixed = T)) %>% 
  # dplyr::mutate(coverage_plot_signal_start = paste0(molecular_trait_id, "&", hit1)) %>% 
  dplyr::mutate(dataset = paste0(dataset, "_", quant)) %>% 
  dplyr::mutate(copy_command = paste0("rsync -r /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_tx_all/",dataset, "/",molecular_trait_id, "\\&", hit1,"* ", gene2,"_idx", idx2, "/", dataset, "_", dplyr::row_number()))

write(sign_coloc_signals_tx$copy_command, file = "copy_tx_coverage_plots.sh")

sign_coloc_signals_txrev = sign_coloc_signals %>% 
  dplyr::filter(quant == "txrev") %>%
  dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = gene1, fixed = T)) %>% 
  # dplyr::mutate(coverage_plot_signal_start = paste0(molecular_trait_id, "&", hit1)) %>% 
  dplyr::mutate(dataset = paste0(dataset, "_", quant)) %>% 
  dplyr::mutate(copy_command = paste0("rsync -r /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_txrev_all/",dataset, "/",molecular_trait_id, "\\&", hit1,"* ", gene2,"_idx", idx2, "/", molecular_trait_id, "_", dplyr::row_number()))

write(sign_coloc_signals_txrev$copy_command, file = "copy_txrev_coverage_plots.sh")

sign_coloc_signals_ge = sign_coloc_signals %>% 
  dplyr::filter(quant == "ge") %>%
  dplyr::mutate(molecular_trait_id = gsub(pattern = ":", replacement = "_", x = gene1, fixed = T)) %>% 
  # dplyr::mutate(coverage_plot_signal_start = paste0(molecular_trait_id, "&", hit1)) %>% 
  dplyr::mutate(dataset = paste0(dataset, "_", quant)) %>% 
  dplyr::mutate(copy_command = paste0("rsync -r /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/",dataset, "/",molecular_trait_id, "___", hit1," ", gene2,"_idx", idx2, "/", molecular_trait_id, "___", hit1,"_", dplyr::row_number()))

write(sign_coloc_signals_ge$copy_command, file = "copy_ge_coverage_plots_new.sh")

sign_coloc_signals$gene2 %>% unique()

library(UpSetR)

# example of list input (list of named vectors)
listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), two = c(1, 2, 4, 5, 
                                                                 10), three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))

# example of expression input
expressionInput <- c(one = 2, two = 1, three = 2, `one&two` = 1, `one&three` = 4, 
                     `two&three` = 1, `one&two&three` = 2)

sign_coloc_signals_summary_upset = sign_coloc_signals_summary_2 %>% 
  dplyr::select(-gene_ids, -gene_names, -count) %>% 
  dplyr::mutate(ge = as.integer(ge)) %>% 
  dplyr::mutate(exon = as.integer(exon)) %>% 
  dplyr::mutate(leafcutter = as.integer(leafcutter)) %>% 
  dplyr::mutate(tx = as.integer(tx)) %>% 
  dplyr::mutate(txrev = as.integer(txrev)) 

list_quants = list()

for (quant in unique(sign_coloc_signals$quant)) {
  list_quants[[quant]] = sign_coloc_signals %>% 
    dplyr::mutate(signal_id = paste0(gene2, "_", idx2)) %>% 
    dplyr::filter(quant == {{quant}}) %>% 
    dplyr::pull(signal_id) %>% 
    unique()
}

upset(fromList(list_quants), order.by = "freq", mainbar.y.label = "Genre Intersections\nID = gene2_idx2")


list_quants = list()

for (quant in unique(sign_coloc_signals$quant)) {
  list_quants[[quant]] = sign_coloc_signals %>% 
    # dplyr::mutate(signal_id = paste0(gene2, "_", idx2)) %>% 
    dplyr::filter(quant == {{quant}}) %>% 
    dplyr::pull(gene2) %>% 
    unique()
}

upset(fromList(list_quants), order.by = "freq", mainbar.y.label = "Genre Intersections\nID = gene2_only")

list_quants = list()
for (quant in unique(sign_coloc_signals$quant)) {
  list_quants[[quant]] = sign_coloc_signals %>% 
    dplyr::mutate(signal_id = paste0(gene2, "_", gene_id)) %>% 
    dplyr::filter(quant == {{quant}}) %>% 
    dplyr::pull(signal_id) %>% 
    unique()
}

upset(fromList(list_quants), order.by = "freq", mainbar.y.label = "Genre Intersections\nID = gene2_geneID")

assigned_colocs = readr::read_tsv("/Users/kerimov/Work/GitHub/update_paper/manuscript_update/coloc/data/Classification_of_colocalising_QTLs.tsv") %>% 
  colSums(na.rm = T) %>% as.data.frame()
assigned_colocs$class = rownames(assigned_colocs)
colnames(assigned_colocs) = c("count", "class")
class_levels = assigned_colocs$class %>% unique()

assigned_colocs = assigned_colocs %>% 
  dplyr::mutate(class = factor(class, levels = class_levels))

# p<-ggplot(data=assigned_colocs, aes(x=class, y=count)) +
#   geom_bar(stat="identity", fill="steelblue")+
#   theme_bw()

all_assigned_signals = readr::read_tsv("/Users/kerimov/Work/GitHub/update_paper/manuscript_update/coloc/data/Classification_of_colocalising_QTLs_all.tsv")[,c(1:17)]

assignments_from_sheets_raw = readr::read_tsv("/Users/kerimov/Work/GitHub/update_paper/manuscript_update/coloc/data/Classification_of_colocalising_QTLs_all_assigned.tsv")
assignments_from_sheets = assignments_from_sheets_raw %>% 
  dplyr::select(gwas_locus = `GWAS locus`, idx2 = `signal (idx2)`, GWAS_variant, gene_count = `gene count`,  eQTL,puQTL,sQTL,apaQTL,ambiguous)

assignments_from_sheets_long = assignments_from_sheets %>% 
  tidyr::pivot_longer(cols = eQTL:ambiguous, names_to = "assigned_to") %>%
  dplyr::filter(!is.na(value)) %>% 
  dplyr::select(-value) %>% 
  dplyr::mutate(assigned_to = factor(assigned_to, levels = class_levels)) %>% 
  dplyr::mutate(single_gene = gene_count == 1)

ggplot(data=assignments_from_sheets_long, aes(x=assigned_to, fill=single_gene)) +
  geom_bar() +
  # geom_text(stat='count', aes(label=..count..))+
  xlab("") +
  ylab("Count of assigned signals") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  theme_classic() +
  theme(legend.position = "none") 

assignments_from_sheets

list_quants_assigned = list()

for (group in unique(assignments_from_sheets_long$assigned_to)) {
  list_quants_assigned[[group]] = assignments_from_sheets_long %>% 
    dplyr::mutate(signal_id = paste0(gwas_locus, "_", idx2)) %>% 
    dplyr::filter(assigned_to == {{group}}) %>% 
    dplyr::pull(signal_id) %>% 
    unique()
}

plot_upset_assigned = upset(fromList(list_quants_assigned), order.by = "freq", mainbar.y.label = "Genre Intersections\nID = gene2_idx2")
ggsave(plot = plot_upset_assigned, filename = "manuscript_update/plots/plot_upset_assigned.pdf", device = "pdf", width = 6, height = 5)

assignments_from_sheets_detected = assignments_from_sheets_raw %>% 
  dplyr::select(gwas_locus = `GWAS locus`, idx2 = `signal (idx2)`, GWAS_variant, unique_gene_id, gene_count = `gene count`,  ge, exon, tx, txrev, leafcutter)

non_ge_colocs_GWAS_loci = assignments_from_sheets_detected %>% 
  dplyr::filter(ge == FALSE) %>% pull(gwas_locus) %>% unique()

  
sign_coloc_signals$quant %>% table

sign_coloc_signals_lc = sign_coloc_signals %>% 
  dplyr::filter(quant == "leafcutter")

sign_coloc_signals_lc$gene2 %>% unique()

sign_coloc_signals %>% 
  dplyr::filter(gene2 %in% non_ge_colocs_GWAS_loci)
  dplyr::filter(quant != "ge") %>% 
  pull(gene1) %>% unique()

# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ ./copy_VidD_sign_colocs_exon.sh
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_exon_all/GTEx_muscle_exon/ENSG00000129933.21_19_19321525_19321903___chr19_19367068_A_C" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_exon_all/TwinsUK_fat_exon/ENSG00000167984.18_16_3539038_3541824___chr16_3620921_C_A" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]
# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ rsync -r '/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_tx_all/Alasoo_2018_macrophage_IFNg_tx/ENST00000576634\&chr16_3539024_AC_A*' 'UKBB.VitD_chr16_3017792-6017591/ENST00000576634_1'
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_tx_all/Alasoo_2018_macrophage_IFNg_tx/ENST00000576634\&chr16_3539024_AC_A*" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]
# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ rsync -r /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_tx_all/Alasoo_2018_macrophage_IFNg_tx/ENST00000576634\&chr16_3539024_AC_A* UKBB.VitD_chr16_3017792-6017591/ENST00000576634_1
# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ chmod +x copy_VidD_sign_colocs_tx.sh
# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ ./copy_VidD_sign_colocs_tx.sh
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_tx_all/Alasoo_2018_macrophage_naive_tx/ENST00000528509&chr11_71459957_G_A*" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_tx_all/GTEx_nerve_tibial_tx/ENST00000430175&chr22_31141547_A_G*" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]
# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ chmod +x copy_VidD_sign_colocs_txrev.sh
# [kerimov@login1] VitD_leafcutter_coloc_coverage_plots_PP4_0.9 $ ./copy_VidD_sign_colocs_txrev.sh
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_txrev_all/GTEx_liver_txrev/ENSG00000145321.grp_2.contained.ENST00000503364&chr4_71742398_C_T*" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]
# rsync: link_stat "/gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_txrev_all/HipSci_iPSC_txrev/ENSG00000070540.grp_1.contained.ENST00000589459&chr17_68452349_A_G*" failed: No such file or directory (2)
# rsync error: some files/attrs were not transferred (see previous errors) (code 23) at main.c(1179) [sender=3.1.2]

### MISSING GE coverage plots
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/BLUEPRINT_SE_monocyte_ge/ENSG00000126214___chr14_103779611_CA_C
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Bossini-Castillo_2019_Treg_naive_ge/ENSG00000172893___chr11_71476248_G_T
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/CommonMind_DLPFC_naive_ge/ENSG00000130518___chr19_18282063_C_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/FUSION_adipose_naive_ge/ENSG00000172890___chr11_71470039_A_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GEUVADIS_LCL_ge/ENSG00000198189___chr4_87341876_AATAG_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_adipose_visceral_ge/ENSG00000132463___chr4_71834529_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_adipose_visceral_ge/ENSG00000134222___chr1_109278889_T_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_adrenal_gland_ge/ENSG00000130518___chr19_18287818_G_C
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_artery_tibial_ge/ENSG00000125629___chr2_118085358_T_C
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_brain_cerebellum_ge/ENSG00000130518___chr19_18275253_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_brain_nucleus_accumbens_ge/ENSG00000130518___chr19_18287818_G_C
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_brain_spinal_cord_ge/ENSG00000130518___chr19_18291332_C_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_colon_sigmoid_ge/ENSG00000198189___chr4_87362944_T_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_colon_transverse_ge/ENSG00000130518___chr19_18277227_C_T
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_esophagus_gej_ge/ENSG00000130518___chr19_18282063_C_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_muscle_ge/ENSG00000125629___chr2_118088243_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_pancreas_ge/ENSG00000134243___chr1_109274570_A_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_skin_not_sun_exposed_ge/ENSG00000108784___chr17_42446094_C_T
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_stomach_ge/ENSG00000134222___chr1_109275684_G_T
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_testis_ge/ENSG00000140678___chr16_30934566_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/GTEx_thyroid_ge/ENSG00000108784___chr17_42445537_A_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/HipSci_iPSC_ge/ENSG00000167394___chr16_30982502_A_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Lepik_2017_blood_ge/ENSG00000143537___chr1_155092683_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Quach_2016_monocyte_naive_ge/ENSG00000164062___chr3_49393727_T_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Schmiedel_2018_CD8_T-cell_anti-CD3-CD28_ge/ENSG00000103351___chr16_3539024_AC_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Schmiedel_2018_CD8_T-cell_naive_ge/ENSG00000103351___chr16_3586808_G_C
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Schmiedel_2018_Th17_memory_ge/ENSG00000160691___chr1_154996879_TCTC_T
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/Schmiedel_2018_Th2_memory_ge/ENSG00000160691___chr1_154972571_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_fat_ge/ENSG00000051108___chr16_56966784_T_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_fat_ge/ENSG00000125629___chr2_118092883_T_G
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_fat_ge/ENSG00000051108___chr16_56966784_T_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_fat_ge/ENSG00000167984___chr16_3620921_C_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_fat_ge/ENSG00000182704___chr11_76779014_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_fat_ge/ENSG00000130522___chr19_18285867_C_CAA
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_skin_ge/ENSG00000139344___chr12_95986028_G_A
# /gpfs/space/projects/eQTLCatalogue/coverage_plots/V6_artemis_ge_all/TwinsUK_skin_ge/ENSG00000163354___chr1_155011887_G_A
