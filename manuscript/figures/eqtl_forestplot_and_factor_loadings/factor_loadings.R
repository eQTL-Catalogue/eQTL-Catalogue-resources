library(tidyverse)

factor_names = tibble(factor = paste0("Factor", c(1:16)),
                      name = c("Universal", "iPSC", "Skin", "Blood",
                               "LCL", "Lymphocyte", "Monocyte & Macrophage", "Adipose",
                               "T cell (anti-CD3/CD28)", "Brain", "Macrophage", "Monocyte",
                               "BLUEPRINT T cell", "ROSMAP Brain", "Muscle", "Neutrophil" ))

mapped_data = "mapping_sn_spMF_K30_a1900_l11100"

# read loading betas and p-values and 
all_betas = read.table(paste0(mapped_data, "_Loadings_beta.txt"))
# to colour significant loadings
pvalues = read.table(paste0(mapped_data, "_Loadings_pvalue_BH.txt"))

loadings_to_tibble <- function(loadings){
  # function to transform sn-spMF output to more convenient
  colnames(loadings) = paste0("Factor", 1:ncol(loadings))
  genes = sapply(strsplit(rownames(loadings), split = " "), "[[", 1)
  variants = sapply(strsplit(rownames(loadings), split = " "), "[[", 2)
  eqtl_ids = paste(variants, genes, sep=".")
  loadings$eqtl_id = eqtl_ids
  loadings$variant = variants
  loadings$gene = genes
  return(dplyr::as_tibble(loadings))
}

get_loadings = function(effects, value_lbl = "Loadings", eqtl = "chr2_160468964_A_T.ENSG00000153250"){
  # transform data
  effects = loadings_to_tibble(effects)
  effect = effects %>% dplyr::filter(eqtl_id == eqtl)
  loadings = effect %>% dplyr::select(starts_with("Factor"))
  # assign names to factors and transform dataframe shape
  colnames(loadings) = factor_names$name
  factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = value_lbl)
}

pvalues = get_loadings(pvalues, "Pvalue")
betas = get_loadings(all_betas)

factors = dplyr::left_join(betas, pvalues)
factors = factors %>% dplyr::mutate(pvalue = ifelse(Pvalue <= 0.05, "<0.05", ">0.05"))

plt = ggplot(factors, aes(x=factor(Factors, levels = Factors), y=Loadings, fill=pvalue)) + 
  geom_col() + 
  xlab("Factors") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.grid = element_blank(), legend.position = 'none')

ggsave("RBMS1_factor_loadings.pdf", plt, width = 3.5, height = 3)

