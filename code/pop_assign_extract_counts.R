library("dplyr")
library("readr")

# NOTE: This script works based on the results of running pop_assign for the VCF files 
#
# EXAMPLE Nextflow script for assign populations of one study
# ==============================================================================
# nextflow run main.nf -profile tartu_hpc -resume\
# -entry pop_assign_only\
# --study_name van_de_Bunt_2015 \
# --vcf_file /gpfs/space/projects/van_de_Bunt_2015/genotypes/genimpute_190921/minimac_out/filtered/van_de_Bunt_2015.filtered.vcf.gz \
# --outdir results/van_de_Bunt_2015
# ==============================================================================

#Import mappings
study_names = list.files(path = "/Users/kerimov/Work/temp_files/pop_assign_7Apr/results")

# Read all the needed sample metadata files 
sample_meta_all <- data.frame()
for (study in study_names) {
  message(study)
  sample_meta <- readr::read_tsv(paste0("/Users/kerimov/Work/GitHub/SampleArcheology/studies/cleaned/",study,".tsv")) %>% 
    dplyr::select(sample_id, genotype_id, rna_qc_passed, genotype_qc_passed, study) %>% 
    dplyr::filter(rna_qc_passed, genotype_qc_passed)

  sample_meta_all <- sample_meta_all %>% 
    rbind(sample_meta)
}

sample_meta_all_qced_genotypes <- sample_meta_all %>% 
  dplyr::pull(genotype_id) %>% 
  unique()

pop_assign_counts_all <- data.frame()
for (study in study_names) {
  message(study)
  pop_assign_result <- readr::read_tsv(paste0("/Users/kerimov/Work/temp_files/pop_assign_7Apr/results/",study,"/",study,"/pop_assign/pop_assigned_abs_0.02_rel_1.7.tsv")) %>% 
    dplyr::mutate(study = study)

  pop_assign_counts_all <- pop_assign_counts_all %>% rbind(pop_assign_result)  
}

pop_assign_counts_all_clean <- pop_assign_counts_all %>% 
  dplyr::filter(genotype_id %in% sample_meta_all_qced_genotypes)

# summarise QCed samples based on relative threshold (1.7) per study
pop_assign_counts_all_clean_counts_rel <- pop_assign_counts_all_clean %>% 
  dplyr::group_by(study, pop_assign_rel_thresh) %>% 
  dplyr::summarise(pop_assign_rel_thresh_count = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(names_from = pop_assign_rel_thresh, values_from = pop_assign_rel_thresh_count) %>% 
  dplyr::select(study, EUR, AFR, SAS, EAS, Admixed) %>% 
  dplyr::rename("Unassigned" = "Admixed")
  
pop_assign_counts_all_clean_counts_rel[is.na(pop_assign_counts_all_clean_counts_rel)] <- 0
pop_assign_counts_all_clean_counts_rel <- pop_assign_counts_all_clean_counts_rel %>% 
  dplyr::mutate(count = rowSums(pop_assign_counts_all_clean_counts_rel %>% select(-study)))

total_samples <- pop_assign_counts_all_clean_counts_rel$count %>% sum()

totals <- as.data.frame(colSums(pop_assign_counts_all_clean_counts_rel %>% select(-study)) %>% t()) %>% 
  dplyr::mutate(study = "Total assigned population across studies")
pop_assign_counts_all_clean_counts_rel <- pop_assign_counts_all_clean_counts_rel %>% 
  rbind(totals) %>% 
  dplyr::mutate(proportion = round(count/total_samples, digits = 4))

write_tsv(pop_assign_counts_all_clean_counts_rel, "../data_tables/pop_assign/population_assignments_rel1.7_thres_QCed.tsv")

##############################################################################################################
# Count the overall assigned populations of unique QC passed samples across all studies (Relative)
##############################################################################################################
# Summarise unique samples across studies (relative)
pop_assign_counts_all_clean_across_studies_rel <- pop_assign_counts_all_clean %>% 
  dplyr::select(-study) %>% 
  dplyr::group_by(genotype_id) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(pop_assign_rel_thresh) %>% 
  dplyr::summarise(pop_assign_rel_thresh_count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(pop_assign_rel_thresh = ifelse(pop_assign_rel_thresh == "Admixed", "Unassigned", pop_assign_rel_thresh)) %>% 
  dplyr::mutate(proportion = round(pop_assign_rel_thresh_count / sum(pop_assign_rel_thresh_count), digits = 3))

total = data.frame(pop_assign_rel_thresh = "Total", 
                   pop_assign_rel_thresh_count = sum(pop_assign_counts_all_clean_across_studies_rel$pop_assign_rel_thresh_count),
                   proportion = sum(pop_assign_counts_all_clean_across_studies_rel$proportion))

pop_assign_counts_all_clean_across_studies_rel <- pop_assign_counts_all_clean_across_studies_rel %>% 
  rbind(total)

readr::write_tsv(pop_assign_counts_all_clean_across_studies_rel, "../data_tables/pop_assign/pop_assign_counts_all_clean_across_studies_rel.tsv")


# summarise QCed samples based on absolute threshold (0.02)
pop_assign_counts_all_clean_counts_abs <- pop_assign_counts_all_clean %>% 
  dplyr::group_by(study, pop_assign_abs_thresh) %>% 
  dplyr::summarise(pop_assign_abs_thresh_count = n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(names_from = pop_assign_abs_thresh, values_from = pop_assign_abs_thresh_count) %>% 
  dplyr::select(study, EUR, AFR, SAS, EAS, Admixed) %>% 
  dplyr::rename("Unassigned" = "Admixed")

pop_assign_counts_all_clean_counts_abs[is.na(pop_assign_counts_all_clean_counts_abs)] <- 0
pop_assign_counts_all_clean_counts_abs <- pop_assign_counts_all_clean_counts_abs %>% 
  dplyr::mutate(count = rowSums(pop_assign_counts_all_clean_counts_abs %>% select(-study)))

total_samples <- pop_assign_counts_all_clean_counts_abs$count %>% sum()

totals <- as.data.frame(colSums(pop_assign_counts_all_clean_counts_abs %>% select(-study)) %>% t()) %>% 
  dplyr::mutate(study = "Total assigned population across studies")
pop_assign_counts_all_clean_counts_abs <- pop_assign_counts_all_clean_counts_abs %>% 
  rbind(totals) %>% 
  dplyr::mutate(proportion = round(count/total_samples, digits = 4))
write_tsv(pop_assign_counts_all_clean_counts_abs, "../data_tables/pop_assign/population_assignments_abs0.02_thres_QCed.tsv")


##############################################################################################################
# Count the overall assigned populations of unique QC passed samples across all studies (Abdolute)
##############################################################################################################

# Summarise unique samples across studies
pop_assign_counts_all_clean_across_studies_abs <- pop_assign_counts_all_clean %>% 
  dplyr::select(-study) %>% 
  dplyr::group_by(genotype_id) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(pop_assign_abs_thresh) %>% 
  dplyr::summarise(pop_assign_abs_thresh_count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(pop_assign_abs_thresh = ifelse(pop_assign_abs_thresh == "Admixed", "Unassigned", pop_assign_abs_thresh)) %>% 
  dplyr::mutate(proportion = round(pop_assign_abs_thresh_count / sum(pop_assign_abs_thresh_count), digits = 3))

total = data.frame(pop_assign_abs_thresh = "Total", 
                   pop_assign_abs_thresh_count = sum(pop_assign_counts_all_clean_across_studies_abs$pop_assign_abs_thresh_count),
                   proportion = sum(pop_assign_counts_all_clean_across_studies_abs$proportion))

pop_assign_counts_all_clean_across_studies_abs <- pop_assign_counts_all_clean_across_studies_abs %>% 
  rbind(total)

readr::write_tsv(pop_assign_counts_all_clean_across_studies_abs, "../data_tables/pop_assign/pop_assign_counts_all_clean_across_studies_abs.tsv")
