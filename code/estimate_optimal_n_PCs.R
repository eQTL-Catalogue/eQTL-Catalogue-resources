library("dplyr")
library("arrow")

read_permuted_pq = function(path){
  data = arrow::read_parquet(path) %>%
    dplyr::mutate(p_fdr = p.adjust(p_beta, method = "fdr"))
}

#Liver, ge, n = 261
liver_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n50/sumstats/QTD000266/QTD000266.permuted.parquet")
liver_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n30/sumstats/QTD000266/QTD000266.permuted.parquet")
liver_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n15/sumstats/QTD000266/QTD000266.permuted.parquet")
liver_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n6/sumstats/QTD000266/QTD000266.permuted.parquet")

dplyr::filter(liver_n50_ge, p_fdr < 0.05)
dplyr::filter(liver_n30_ge, p_fdr < 0.05)
dplyr::filter(liver_n15_ge, p_fdr < 0.05)
dplyr::filter(liver_n6_ge, p_fdr < 0.05)
#Optimal: n = 50, n = 30 not far behind

#Liver, lc, n = 261
liver_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n50/sumstats/QTD000270/QTD000270.permuted.parquet")
liver_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n30/sumstats/QTD000270/QTD000270.permuted.parquet")
liver_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n15/sumstats/QTD000270/QTD000270.permuted.parquet")
liver_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_liver_v10_ge_leaf_n6/sumstats/QTD000270/QTD000270.permuted.parquet")

dplyr::filter(liver_n50_ge, p_fdr < 0.05)
dplyr::filter(liver_n30_ge, p_fdr < 0.05)
dplyr::filter(liver_n15_ge, p_fdr < 0.05)
dplyr::filter(liver_n6_ge, p_fdr < 0.05)
#Optimal: n = 6, but n = 15 and n = 30 are similar

#Blood, ge, n = 800
blood_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n50/sumstats/QTD000356/QTD000356.permuted.parquet")
blood_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n30/sumstats/QTD000356/QTD000356.permuted.parquet")
blood_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n15/sumstats/QTD000356/QTD000356.permuted.parquet")
blood_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n6/sumstats/QTD000356/QTD000356.permuted.parquet")

dplyr::filter(blood_n50_ge, p_fdr < 0.05)
dplyr::filter(blood_n30_ge, p_fdr < 0.05)
dplyr::filter(blood_n15_ge, p_fdr < 0.05)
dplyr::filter(blood_n6_ge, p_fdr < 0.05)
#Optimal, n = 50, n = 30 significantly worse

#Blood, lc, n = 800
blood_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n50/sumstats/QTD000360/QTD000360.permuted.parquet")
blood_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n30/sumstats/QTD000360/QTD000360.permuted.parquet")
blood_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n15/sumstats/QTD000360/QTD000360.permuted.parquet")
blood_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_blood_v10_ge_leaf_n6/sumstats/QTD000360/QTD000360.permuted.parquet")

dplyr::filter(blood_n50_ge, p_fdr < 0.05)
dplyr::filter(blood_n30_ge, p_fdr < 0.05)
dplyr::filter(blood_n15_ge, p_fdr < 0.05)
dplyr::filter(blood_n6_ge, p_fdr < 0.05)
#Optimal, n = 50, but all others are very similar including n = 6.

#Kidney, ge, n = 103
kidney_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n50/sumstats/QTD000261/QTD000261.permuted.parquet")
kidney_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n30/sumstats/QTD000261/QTD000261.permuted.parquet")
kidney_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n15/sumstats/QTD000261/QTD000261.permuted.parquet")
kidney_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n6/sumstats/QTD000261/QTD000261.permuted.parquet")

dplyr::filter(kidney_n50_ge, p_fdr < 0.05)
dplyr::filter(kidney_n30_ge, p_fdr < 0.05)
dplyr::filter(kidney_n15_ge, p_fdr < 0.05)
dplyr::filter(kidney_n6_ge, p_fdr < 0.05)
#Optimal, n = 15, n = 30 worse but not too much

#Kidney, leafcutter, n = 103
kidney_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n50/sumstats/QTD000265/QTD000265.permuted.parquet")
kidney_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n30/sumstats/QTD000265/QTD000265.permuted.parquet")
kidney_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n15/sumstats/QTD000265/QTD000265.permuted.parquet")
kidney_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_kidney_v10_ge_leaf_n6/sumstats/QTD000265/QTD000265.permuted.parquet")

dplyr::filter(kidney_n50_ge, p_fdr < 0.05)
dplyr::filter(kidney_n30_ge, p_fdr < 0.05)
dplyr::filter(kidney_n15_ge, p_fdr < 0.05)
dplyr::filter(kidney_n6_ge, p_fdr < 0.05)
#Optimal, n = 6, n = 15 worse but not too much

#Breast, ge, n = 511
breast_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n50/sumstats/QTD000211/QTD000211.permuted.parquet")
breast_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n30/sumstats/QTD000211/QTD000211.permuted.parquet")
breast_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n15/sumstats/QTD000211/QTD000211.permuted.parquet")
breast_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n6/sumstats/QTD000211/QTD000211.permuted.parquet")

dplyr::filter(breast_n50_ge, p_fdr < 0.05)
dplyr::filter(breast_n30_ge, p_fdr < 0.05)
dplyr::filter(breast_n15_ge, p_fdr < 0.05)
dplyr::filter(breast_n6_ge, p_fdr < 0.05)
#Optimal, n = 6, n = 15 worse but not too much

#Breast, ge, n = 511
breast_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n50/sumstats/QTD000211/QTD000211.permuted.parquet")
breast_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n30/sumstats/QTD000211/QTD000211.permuted.parquet")
breast_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n15/sumstats/QTD000211/QTD000211.permuted.parquet")
breast_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n6/sumstats/QTD000211/QTD000211.permuted.parquet")

dplyr::filter(breast_n50_ge, p_fdr < 0.05)
dplyr::filter(breast_n30_ge, p_fdr < 0.05)
dplyr::filter(breast_n15_ge, p_fdr < 0.05)
dplyr::filter(breast_n6_ge, p_fdr < 0.05)
#Optimal, n = 50, n = 30 is worse

#Breast, lc, n = 511
breast_n50_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n50/sumstats/QTD000215/QTD000215.permuted.parquet")
breast_n30_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n30/sumstats/QTD000215/QTD000215.permuted.parquet")
breast_n15_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n15/sumstats/QTD000215/QTD000215.permuted.parquet")
breast_n6_ge = read_permuted_pq("large_data/GTEx_n_pheno_pcs/GTEx_breast_v10_ge_leaf_n6/sumstats/QTD000215/QTD000215.permuted.parquet")

dplyr::filter(breast_n50_ge, p_fdr < 0.05)
dplyr::filter(breast_n30_ge, p_fdr < 0.05)
dplyr::filter(breast_n15_ge, p_fdr < 0.05)
dplyr::filter(breast_n6_ge, p_fdr < 0.05)
#Optimal, n = 30


#Final recommendations
#If n <= 150, ge, exon n_pheno_pcs = 15, tx, txrev,lc trait n_pheno_pcs = 6
#If 150 < n <= 350, ge, exon n_pheno_pcs = 30, tx, txrev,lc n_pheno_pcs = 15
#if n > 350, ge, exon n_pheno_pcs = 50, tx, txrev,lc n_pheno_pcs = 30


