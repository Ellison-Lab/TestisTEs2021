library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)
library(ggtext)

source("workflow/fig-scripts/theme.R")

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr/", format='arrow')

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+"))))

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
  filter(qval <  optimal_ica[['qval']]) %>%
  collect() %>%
  filter(as.character(module)==tep.name) %>%
  pull(X1)

tep.tes <- tep.members[!str_detect(tep.members,'FBgn')]

tep_te_expr_df <- map_df(as.list(as.numeric(unique(pull(collect(w1118.obs),'clusters')))) %>% set_names(.,.),
                         ~{filter(w1118.expr, (clusters == .) & (gene_id %in% tep.tes)) %>% collect()}) %>%
  dplyr::select(index, gene_id, expression, clusters) %>%
  mutate(clusters = as.character(clusters)) %>%
  mutate(umis = exp(expression) - 1)

tep_gene_expr_df <- map_df(as.list(as.numeric(unique(pull(collect(w1118.obs),'clusters')))) %>% set_names(.,.),
                         ~{filter(w1118.expr, (clusters == .) & (gene_id %in% c("FBgn0036470"))) %>% collect()}) %>%
  dplyr::select(index, gene_id, expression, clusters) %>%
  mutate(clusters = as.character(clusters)) %>%
  mutate(umis = exp(expression) - 1)

top_tep_te_rnk_df <- tep_te_expr_df %>%
  left_join(collect(w1118.obs), by=c(index='X1',clusters='clusters')) %>%
  group_by(gene_id, clusters2, clusters) %>%
  summarize(mean_expr = mean(expression)) %>%
  ungroup() %>%
  group_by(clusters2) %>%
  mutate(rnk = dense_rank(-mean_expr)) %>%
  arrange(clusters2,rnk) %>%
  filter(rnk <= 20) %>%
  left_join(rename.table)

top.expressers <- top_tep_te_rnk_df %>% group_by(clusters.rename) %>%
  summarise(mean_expr = mean(mean_expr)) %>%
  arrange(-mean_expr) %>%
  pull(clusters.rename) %>%
  head(1)

df <- tep_te_expr_df %>%
  left_join(collect(w1118.obs), by=c(index='X1',clusters='clusters')) %>%
  left_join(rename.table) %>%
  filter(clusters.rename == top.expressers) %>%
  filter(gene_id %in% c(top_tep_te_rnk_df$gene_id,tep_gene_expr_df$gene_id))
  
g0 <- df %>% filter(!str_detect(gene_id,"FBgn")) %>%
  ggplot( aes(reorder(gene_id,expression),expression)) +
  geom_boxplot(fill="lightgray") +
  theme_gte21() +
  facet_wrap(~clusters.rename,scales='free') +
  ylab('log-norm UMIs') + xlab('') +
  theme(aspect.ratio = 2) +
  guides(fill=F) +
  coord_flip()

g <- tep_te_expr_df %>%
  bind_rows(tep_gene_expr_df) %>%
  left_join(collect(w1118.obs), by=c(index='X1',clusters='clusters')) %>%
  left_join(rename.table) %>%
  filter(gene_id %in% c("QUASIMODO2","ACCORD2","FBgn0036470")) %>%
  mutate(gene_id = ifelse(gene_id == "FBgn0036470", "EAChm",gene_id)) %>%
  mutate(gene_id=fct_relevel(gene_id, "QUASIMODO2","ACCORD2","EAChm")) %>%
  ggplot(aes(clusters.rename, expression, fill=clusters.rename)) +
  geom_violin(scale = "width") +
  scale_fill_gte21() +
  theme_gte21()  +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  theme(legend.text = element_text(size=7/.pt), legend.title = element_text(size=7/.pt), strip.text.x = element_markdown()) +
  xlab('') + ylab('') +
  ylab('log-norm UMIs') + xlab('') +guides(fill=F) +
  facet_wrap(~gene_id, ncol=1)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
