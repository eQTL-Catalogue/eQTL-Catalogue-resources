Nathan_2022 = readr::read_table("big_data/popAssign_tables/Nathan_2022_pop_assigned_abs_0.02_rel_1.7.tsv")
MacroMap = readr::read_table("big_data/popAssign_tables/MacroMap_pop_assigned_abs_0.02_rel_1.7.tsv")
Cytoimmgen = readr::read_table("big_data/popAssign_tables/Cytoimmgen_pop_assigned_abs_0.02_rel_1.7.tsv")
Perez1 = readr::read_table("big_data/popAssign_tables/Perex_2022_array1_pop_assigned_abs_0.02_rel_1.7.tsv")
Perez2 = readr::read_table("big_data/popAssign_tables/Perex_2022_array2_pop_assigned_abs_0.02_rel_1.7.tsv")
OneK1K = readr::read_table("big_data/popAssign_tables/OneK1K_pop_assigned_abs_0.02_rel_1.7.tsv")

GAINS = readr::read_table("big_data/popAssign_tables/GAINS_pop_assigned_abs_0.02_rel_1.7.tsv")
d = readr::read_table("~/projects/SampleArcheology/studies/cleaned/GAinS.tsv") %>%
  dplyr::filter(genotype_qc_passed, rna_qc_passed)
gains2 = dplyr::filter(GAINS, genotype_id %in% d$genotype_id)

AFR_LCL = readr::read_table("big_data/popAssign_tables/AFR_LCL_pop_assigned_abs_0.02_rel_1.7.tsv")

PISA = readr::read_tsv("big_data/popAssign_tables/PISA_pop_assigned_abs_0.02_rel_1.7.tsv")

Nassiri_2025 = readr::read_tsv("big_data/popAssign_tables/Nassiri_2025_pop_assigned_abs_0.02_rel_1.7.tsv")

d = readr::read_tsv("~/projects/SampleArcheology/studies/cleaned/MAGE.tsv")

Aygun_2021 = readr::read_tsv("big_data/popAssign_tables/Aygun_2021_pop_assigned_abs_0.02_rel_1.7.tsv")

Walker_2019 = readr::read_tsv("big_data/popAssign_tables/Walker_2019_pop_assigned_abs_0.02_rel_1.7.tsv")

IBDverse = readr::read_tsv("big_data/popAssign_tables/IBDverse_pop_assigned_abs_0.02_rel_1.7.tsv")

DICE = readr::read_tsv("big_data/popAssign_tables/DICE_popassign.tsv")

GTEx_v10 = readr::read_tsv("big_data/popAssign_tables/GTEx_pop_assigned_abs_0.02_rel_1.7.tsv")

INTERVAL = readr::read_tsv("big_data/popAssign_tables/INTERVAL_pop_assigned_abs_0.02_rel_1.7.tsv")

Natri_2024 = readr::read_tsv("big_data/popAssign_tables/Natri_2024_pop_assigned_abs_0.02_rel_1.7.tsv")

Sun_2018 = readr::read_tsv("big_data/popAssign_tables/Sun_2018_pop_assigned_abs_0.02_rel_1.7.tsv")
