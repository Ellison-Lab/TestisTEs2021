library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

te.lookup <- read_tsv('~/work/TestisTpn/data/te_id_lookup.curated.tsv.txt')


w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.gep_usage <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_usage/", format='arrow')
w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

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

tep.members <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  filter(as.character(module)==tep.name) %>%
  pull(X1)

umap_labs <- w1118.obs %>%
  group_by(clusters2) %>%
  collect() %>%
  summarise(x = mean(X_umap1), y=mean(X_umap2)) %>%
  ungroup() %>%
  left_join(rename.table)

df <- w1118.gep_usage %>%
  filter(cons == tep.name) %>%
  collect() %>%
  left_join(collect(w1118.obs)) %>%
  arrange(usage)

umap_labs <- w1118.obs %>%
  group_by(clusters2) %>%
  collect() %>%
  summarise(x = mean(X_umap1), y=mean(X_umap2)) %>%
  ungroup() %>%
  left_join(rename.table)

g <-  ggplot(df, aes(X_umap1,X_umap2)) +
  geom_point(size=0.75, aes(color=usage)) +
  scale_color_gradient2(low='steelblue', mid='lightgray',high='tomato', midpoint = 0.04) +
  theme_classic() +
  coord_fixed()  +
  geom_text(data=umap_labs, aes(x + sign(5-x) *2, y, label=clusters.rename), face='bold', size=rel(3)) +
  theme(aspect.ratio = 1) +
  xlab("UMAP1") + ylab("UMAP2")

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])

