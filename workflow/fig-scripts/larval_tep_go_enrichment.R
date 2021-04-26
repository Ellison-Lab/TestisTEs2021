library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)

source("workflow/fig-scripts/theme.R")

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format='arrow')
w1118.gep_enr <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_enr", format='arrow')
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

df <- w1118.gep_enr %>%
  collect() %>%
  filter(cluster == tep.name) %>%
  group_by(ont) %>%
  slice_max(score, n=10) %>%
  mutate(rnk = dense_rank(score))

g <-  df %>%
  ggplot(aes(score, reorder(Term, rnk))) +
  geom_col(aes(fill=ont)) +
  facet_wrap(~ont, scales = "free", ncol = 1) +
  theme_gte21() +
  scale_fill_brewer(type='qual', palette = 3) +
  theme(aspect.ratio = 1) +
  xlab("-log10(pval)") + ylab("") +
  guides(fill=F)


agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()


saveRDS(g,snakemake@output[['ggp']])

write_tsv(df,snakemake@output[['dat']])

