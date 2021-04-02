library(tidyverse)
library(arrow)
library(ragg)
library(ggpubr)

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow') %>% collect()
w1118.te_umis <- open_dataset("results/finalized/larval-w1118-testes/umis/", format='arrow') %>% collect()

df <- w1118.obs %>% dplyr::select(X1, clusters2) %>%
  left_join(w1118.te_umis, .,by=c(cell="X1")) %>%
  left_join(rename.table)  

df.per_cell <- df %>% group_by(batch, clusters.rename, cell) %>%
  summarize(UMIs=sum(UMIs), .groups = "drop")

g <- ggplot(df.per_cell, aes(clusters.rename,UMIs)) +
  #geom_boxplot() +
  geom_violin(aes(fill=clusters.rename),draw_quantiles = c(0.5),scale = 'width') +
  theme_classic() +
  xlab("") + ylab('UMIs') +
  guides(fill=F) +
  scale_fill_brewer(type='qual', palette = 8) +
  theme(aspect.ratio = 0.3) +
  theme(plot.caption= element_text(hjust=0.5, face='italic', size=rel(1.2)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text.y = element_text(size=rel(1)),
        axis.text.x = element_text(size=rel(1.5), angle=90, hjust=1, vjust=0.5)) +
  stat_compare_means()

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
