library(tidyverse)
library(ggplot2)

sizes = readr::read_tsv("Supplementary_Tables_R3 - pop_assign_summary.tsv")
sizes = sizes %>% filter(`Assigned population` != "TOTAL")
sizes = sizes %>% mutate(`Assigned population` = gsub("Admixed", "NA", `Assigned population`))
sizes$`Assigned population` = factor(sizes$`Assigned population`, levels=sizes$`Assigned population`)

plt = ggplot(sizes, aes(x=`Assigned population`, `Sample Size`))+
  geom_col() +
  xlab("Assigned population") + 
  geom_text(aes(label = Percent), vjust = -0.5) + 
  theme_light() + 
  theme(panel.grid = element_blank())

ggsave("assigned_populations.pdf", width=5, height=5)