library(tidyverse)
library(arrow)
library(ragg)
library(ggrepel)

source("workflow/fig-scripts/theme.R")

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+"))))

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

umap_labs <- w1118.obs %>%
  group_by(clusters2) %>%
  collect() %>%
  summarise(x = mean(X_umap1), y=mean(X_umap2)) %>%
  ungroup() %>%
  left_join(rename.table)

g <- collect(w1118.obs) %>%
  left_join(rename.table) %>%
  ggplot(aes(X_umap1, X_umap2)) +
  geom_point(aes(color=reorder(clusters.rename, clusters.rename)), size=rel(0.75)) +
  scale_color_gte21() +
  theme_gte21() +
  labs(color="") +
  geom_text_repel(data=umap_labs, aes(x + sign(5-x) *2, y, label=clusters.rename), face='bold', size=rel(3)) +
  xlab("UMAP1") + ylab("UMAP2") +
  theme(plot.caption= element_text(hjust=0.5, face='italic', size=rel(1.2)), axis.title = element_text(size = rel(1.2)), axis.text = element_text(size=rel(1.5))) +
  coord_fixed() +
  guides(color=F)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(collect(w1118.obs),snakemake@output[['dat']])
