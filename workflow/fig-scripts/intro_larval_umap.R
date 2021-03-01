library(tidyverse)
library(arrow)
library(ragg)

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

umap_labs <- w1118.obs %>%
  group_by(clusters2) %>%
  collect() %>%
  summarise(x = mean(X_umap1), y=mean(X_umap2)) %>%
  ungroup() %>%
  left_join(rename.table)


g <- collect(w1118.obs) %>%
  ggplot(aes(X_umap1, X_umap2)) +
  geom_point(aes(color=clusters2), size=rel(0.75)) +
  scale_color_brewer(type='qual',palette='Set3', name='') +
  theme_classic() +
  #coord_fixed(xlim = c(-15, 8)) +
  theme(legend.position = 'right') +
  labs(color="") +
  geom_text(data=umap_labs, aes(x + sign(x) *1, y, label=clusters.rename), face='bold', size=rel(3)) +
  theme(plot.caption= element_text(hjust=0.5, face='italic')) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(plot.caption= element_text(hjust=0.5, face='italic', size=rel(1.2)), axis.title = element_text(size = rel(1.2)), axis.text = element_text(size=rel(1.5)))

#ggsave(filename = snakemake@output[['png']], plot=g, height=10, width = 10, scale=2)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(collect(w1118.obs),snakemake@output[['dat']])
