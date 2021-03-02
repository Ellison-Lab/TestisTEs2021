library(tidyverse)
library(arrow)
library(ragg)

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr", format='arrow')

df <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.list %>% set_names(.,.),
       ~{filter(w1118.expr, clusters == . & gene_id %in% unique(te.lookup$merged_te)) %>% collect()}) %>%
  dplyr::select(index, gene_id, expression) %>%
  mutate(expression = exp(expression) - 1) %>%
  group_by(index) %>%
  summarize(expression=log2(sum(expression) + 1)) %>%
  ungroup() %>%
  left_join(collect(w1118.obs), by=c(index='X1'))
  
df <- df %>%
  left_join(rename.table)  


g <- ggplot(df, aes(clusters.rename,expression)) +
  geom_violin(aes(fill=clusters.rename),draw_quantiles = c(0.5),scale = 'width') +
  theme_classic() +
  xlab("") + ylab('log2(Normalized TE UMIs + 1)') +
  guides(fill=F) +
  scale_fill_brewer(type='qual', palette = 8) +
  theme(aspect.ratio = 0.3) +
  theme(plot.caption= element_text(hjust=0.5, face='italic', size=rel(1.2)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text.y = element_text(size=rel(1)),
        axis.text.x = element_text(size=rel(1.5), angle=90, hjust=1, vjust=0.5))

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
