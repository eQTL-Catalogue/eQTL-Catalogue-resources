library(tidyverse)
library(ggplot2)

sizes = readr::read_tsv("Supplementary_Tables_R3 - pop_assign_summary.tsv")
sizes = sizes %>% filter(`Assigned population` != "TOTAL")
sizes = sizes %>% mutate(`Assigned population` = gsub("Admixed", "NA", `Assigned population`))
sizes$`Assigned population` = factor(sizes$`Assigned population`, levels=sizes$`Assigned population`)

plt = ggplot(sizes, aes(x=`Assigned population`, `Sample Size`))+
  geom_col() +
  xlab("Assigned population") + 
  ylab("Sample size") +
  geom_text(aes(label = paste0(Percent,"%")), vjust = -0.5, size = 3) + 
  theme_light() + 
  theme(panel.grid = element_blank()) +
  ylim(0,5700)

ggsave("assigned_populations.pdf", width=2.8, height=2.8)
