library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)

source("workflow/fig-scripts/theme.R")

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

tep.name <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  dplyr::select(module, X1, repClass, repFamily) %>%
  distinct() %>%
  group_by(module) %>%
  summarize(n_tes=n()) %>%
  arrange(-n_tes) %>%
  head(1) %>%
  pull(module) %>% as.character()

df <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  group_by(module) %>% tally()

g <- ggplot(df, aes(n)) +
  geom_histogram(aes(fill=..count..), color='black') +
  theme_gte21() +
  scale_y_continuous(expand = expand_scale(mult = c(0,0), add = c(0,1))) +
  scale_fill_fermenter(direction = 1, palette = 13) +
  scale_fill_gte21("diverging.wide", discrete = F, reverse = T) +
  theme(aspect.ratio = 0.5, plot.caption = element_text(hjust=0)) +
  xlab('module sizes')

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
