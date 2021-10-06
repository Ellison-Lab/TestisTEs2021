library(tidyverse)
library(GenomicRanges)
library(arrow)
library(jsonlite)
library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)

source("workflow/fig-scripts/theme.R")

# -----------------------------------------------------------
# reviewer 1 count request
# -----------------------------------------------------------

te.lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

tes <- te.lookup$gene_id %>% unique

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

geps <- open_dataset('results/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

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

top_mods <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
  arrange(-n_tes)

tep_tes <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  filter(module == tep.name) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  pull(X1) %>%
  unique()

lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")


te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr", format='arrow')

df <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.numeric %>% as.list %>% set_names(.,.),
             ~{filter(w1118.expr, clusters == . & gene_id %in% unique(te.lookup$merged_te)) %>% collect()}) %>%
  dplyr::select(index, gene_id, expression) %>%
  mutate(expression = exp(expression) - 1) %>%
  left_join(collect(w1118.obs), by=c(index='X1')) %>%
  left_join(rename.table) %>%
  group_by(clusters.rename,gene_id) %>%
  filter(sum(expression > 0) >= n()/2) %>% # TEs expressed in at least a half of the cluster
  summarise(expression = sum(expression)) %>%
  ungroup()

g <- df %>%
  ggplot(aes(clusters.rename)) +
  geom_bar() +
  theme_gte21() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust = 1)) +
  xlab("scRNA-seq clusters") + ylab("N TEs")

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])



# df %>% 
#   filter(str_detect(clusters.rename,"Pig") | str_detect(clusters.rename,"Termin")) %>%
#   mutate(clusters.rename=droplevels(clusters.rename)) %>% 
#   split(.,.$clusters.rename) %>% 
#   map(~pull(.,gene_id)) %>%
#   {.[["9/Pigment"]][.[["9/Pigment"]] %in% .[["10/TerminalEpithelial"]]]}
# 
# df %>% 
#   filter(str_detect(clusters.rename,"1/Spermatogonia")) %>%
#   pull(gene_id)
# 
# df %>% 
#   filter(str_detect(clusters.rename,"Pigment")) %>%
#   pull(gene_id)
# 
# df %>% 
#   filter(str_detect(clusters.rename,"Terminal")) %>%
#   pull(gene_id)
# 
# df %>% 
#   filter(str_detect(clusters.rename,"Cyst")) %>%
#   pull(gene_id) %>%
#   unique()
# 
# df %>% 
#   filter(str_detect(clusters.rename,"3/Spermatocyte")) %>%
#   pull(gene_id) %>%
#   unique()
# 
# 
# df %>% 
#   filter(str_detect(clusters.rename,"2/TransitionalSpermatocyte")) %>%
#   pull(gene_id) %>%
#   unique()
# 
# df %>% 
#   filter(str_detect(clusters.rename,"2/Trans") | str_detect(clusters.rename,"3/Sperm")) %>%
#   mutate(clusters.rename=droplevels(clusters.rename)) %>% 
#   split(.,.$clusters.rename) %>% 
#   map(~pull(.,gene_id)) %>%
#   {.[["3/Spermatocyte"]][.[["3/Spermatocyte"]] %in% .[["2/TransitionalSpermatocyte"]]]}
