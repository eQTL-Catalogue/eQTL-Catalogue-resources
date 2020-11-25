library("dplyr")
library("data.table")
library("ggplot2")

colocs = readr::read_tsv("merged_colocs.tsv") %>%
  dplyr::filter(qtl_subset != "qtl_subset") %>%
  dplyr::mutate(PP.H4.abf = as.numeric(PP.H4.abf),
                PP.H3.abf = as.numeric(PP.H3.abf)) %>%
  dplyr::mutate(quant = ifelse(qtl_subset %like% "_tx", "transcript", "gene")) %>%
  dplyr::mutate(quant = ifelse(qtl_subset %like% "_txrev", "txrevise", quant)) %>%
  dplyr::mutate(quant = ifelse(qtl_subset %like% "_exon", "exon", quant))


plt = ggplot(colocs, aes(x = PP.H3.abf, y = PP.H4.abf, color = quant)) + 
  geom_point(alpha = 0.8) +
  theme_light() +
  theme(panel.grid = element_blank()) +
  xlab("PP3 (distinct causal variants)") +
  ylab("PP4 (shared causal variant)")
  
ggsave("HMGCR_coloc_plot.pdf", plot = plt, width = 3.5, height = 2.7)




#Make effect size plots
exon_effect_sizes = readr::read_tsv("HMGCR_exon.tsv")

# calculate 95% confidence interval
conf.level = 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
exon_effect_sizes = dplyr::mutate(exon_effect_sizes, interval = ci.value * se)

plt = ggplot(exon_effect_sizes, aes(x = file_name, y = beta, ymin = beta - interval, ymax = beta + interval)) + 
  geom_point() + 
  geom_hline(yintercept=0, colour="grey") + 
  geom_errorbar(width = 0.1) + 
  xlab("Dataset") + 
  ylab("Effect size") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(1, 1, 1, 3, "cm"),panel.grid = element_blank(),)+
  geom_hline(yintercept = 0)

ggsave("HMGCR_exon_forest.pdf", plot = plt, width = 8, height = 3)


#Make gene expression effect size plots
gene_effect_sizes = readr::read_tsv("HMGCR_gene.tsv")

# calculate 95% confidence interval
conf.level = 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
gene_effect_sizes = dplyr::mutate(gene_effect_sizes, interval = ci.value * se)

plt = ggplot(gene_effect_sizes, aes(x = file_name, y = beta, ymin = beta - interval, ymax = beta + interval)) + 
  geom_point() + 
  geom_hline(yintercept=0, colour="grey") + 
  geom_errorbar(width = 0.1) + 
  xlab("Dataset") + 
  ylab("Effect size") +
  theme_light() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(1, 1, 1, 3, "cm"),panel.grid = element_blank(),)+
  geom_hline(yintercept = 0)

ggsave("HMGCR_gene_forest.pdf", plot = plt, width = 8, height = 3)


