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

