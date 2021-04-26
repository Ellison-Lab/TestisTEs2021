library(tidyverse)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr", format='arrow')

markers <- c("piwi",
             "qin",
             "vas",
             "AGO3",
             "armi",
             "aub",
             "csul",
             "Gasz",
             "Hen1",
             "krimp",
             "mino",
             "Nbr",
             "shu",
             "spn-E",
             "vls",
             "vret",
             "zuc",
             "cuff",
             "del",
             "rhi",
             "moon",
             "arx",
             "Panx",
             "nxf2",
             "egg",
             "Su(var)3-3",
             "wde",
             "mael",
             "SoYb",
             "BoYb",
             "tej")

df <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.list %>% set_names(.,.),
             ~{filter(w1118.expr, clusters == . & gene_symbol %in% markers) %>% collect()}) %>%
  dplyr::select(index, gene_symbol, expression) %>%
  mutate(expression = exp(expression) - 1) %>%
  left_join(collect(w1118.obs), by=c(index='X1')) %>%
  distinct()
  #filter(str_detect(clusters2, 'Sperm'))

df <- df %>%
  left_join(rename.table) %>%
  group_by(clusters.rename, gene_symbol) %>%
  summarize(pct.expressing = sum(expression > 0)/n(), mean.expression=mean(expression))

g <- df %>%
  ggplot(aes(gene_symbol, clusters.rename)) +
  geom_point(aes(size=pct.expressing, fill=mean.expression), shape=21) +
  scale_fill_fermenter(palette = 8, direction = 1, name='mean Log1p normalized UMIs', guide=guide_legend(label.position = 'bottom', title.position = 'top')) +
  scale_size(range=c(0, rel(7)), name='Proportion expressing', guide=guide_legend(label.position = 'bottom', title.position = 'top')) +
  theme_gte21()  +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  theme(legend.text = element_text(size=rel(0.5)), legend.title = element_text(size=rel(0.5))) +
  theme(aspect.ratio = 2) +
  coord_flip() +
  xlab('') + ylab('')

agg_png(snakemake@output[['png']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
