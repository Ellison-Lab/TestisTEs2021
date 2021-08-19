library(tidyverse)
library(arrow)
library(ragg)
library(ggpubr)

source("workflow/fig-scripts/theme.R")

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow') %>% collect()

w1118.te_umis <- open_dataset("results/finalized/larval-w1118-testes/umis/", format='arrow') %>% collect()

df <- w1118.obs %>% dplyr::select(X1, clusters2) %>%
  left_join(w1118.te_umis, .,by=c(cell="X1"))

df.per_cell <- df %>% group_by(batch, clusters2, cell) %>%
  summarize(UMIs=sum(UMIs), .groups = "drop")

df.per_cell <- df.per_cell %>%
  left_join(rename.table)  

g <- ggplot(df.per_cell, aes(clusters.rename,UMIs)) +
  geom_violin(aes(fill=clusters.rename),draw_quantiles = c(0.5),scale = 'width') +
  theme_gte21() +
  xlab("") + ylab('TE-derived UMIs (raw)') +
  guides(fill=F) +
  scale_fill_gte21() +
  theme(plot.caption= element_text(hjust=0.5, face='italic'),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.title = element_text(size = rel(1))) +
  stat_compare_means(label.y.npc = 0.9, label.x.npc = 0.1)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

raw.stats <- kruskal.test(y~x,data=g$layers[[2]]$compute_aesthetics(data=g$data, plot=g)) %>%
  broom::tidy()

stats.export <- raw.stats %>%
  mutate(script= "larval_raw_te_umis.R") %>%
  mutate(desc = "raw UMI counts") %>%
  mutate(func = "stats::kruskal.test/ggpubr::stat_compare_means") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "Inter-cluster expression") %>%
  mutate(alternative=NA) %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
