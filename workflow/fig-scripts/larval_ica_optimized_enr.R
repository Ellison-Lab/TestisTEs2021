library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)

source("workflow/fig-scripts/theme.R")

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format='arrow')
w1118.gep_enr <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_enr", format='arrow')

df <- w1118.gep_enr %>%
  collect() %>%
  filter(weight01 < 0.05) %>%
  group_by(ont,Term) %>%
  mutate(is.unique = n()) %>%
  group_by(ont, cluster) %>%
  summarise(has.unique = any(is.unique==1),.groups = "drop_last") %>%
  summarise(pct.unique = sum(has.unique)/n(), n.unique = sum(has.unique))

g <- ggplot(df,aes(ont,pct.unique)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent modules w/ unique enrichent") +
  theme_gte21() +
  theme(axis.title.y = element_text(size=7/.pt), axis.title.x = element_blank())

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()


saveRDS(g,snakemake@output[['ggp']])

write_tsv(df,snakemake@output[['dat']])
