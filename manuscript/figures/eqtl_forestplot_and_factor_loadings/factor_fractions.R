library(tidyverse)

factor_names = tibble(factor = paste0("Factor", c(1:16)),
                      name = c("Universal", "iPSC", "Skin", "Blood",
                               "LCL", "Lymphocyte", "Monocyte & Macrophage", "Adipose",
                               "T cell (anti-CD3/CD28)", "Brain", "Macrophage", "Monocyte",
                               "BLUEPRINT T cell", "ROSMAP Brain", "Muscle", "Neutrophil"))

loadings_file = "sn_spMF_K30_a1900_l11100_Loadings_beta_alpha0.05_corrected.txt.gz"
loadings = read.delim(loadings_file)

loadings_to_tibble <- function(loadings){
  colnames(loadings) = paste0("Factor", 1:ncol(loadings))
  genes = sapply(strsplit(rownames(loadings), split = " "), "[[", 1)
  variants = sapply(strsplit(rownames(loadings), split = " "), "[[", 2)
  eqtl_ids = paste(variants, genes, sep=".")
  loadings$eqtl_id = eqtl_ids
  loadings$variant = variants
  loadings$gene = genes
  return(dplyr::as_tibble(loadings))
}

loadings = loadings_to_tibble(loadings)

n_eqtls = nrow(loadings)

loadings = loadings %>% select(starts_with("Factor"))
colnames(loadings) = factor_names$name

factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = "Loadings")
fractions = factors %>% group_by(Factors) %>% summarise(count = sum(Loadings != 0), Fraction = count / n_eqtls)

plt = ggplot(fractions, aes(reorder(Factors, -Fraction), Fraction)) +
  geom_col() + 
  theme_light() +
  xlab("Factors") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), text = element_text(size=16), panel.grid = element_blank())
ggsave("factor_assignment.png", plot = plt, width = 8.2, height = 4)
