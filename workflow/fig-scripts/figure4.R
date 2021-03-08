library(tidyverse)

library(patchwork)

y_genes <- read_rds('results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds')

late_sperm <- read_rds("results/figs/larval_later_sperm_marker_umis/larval_later_sperm_marker_umis.ggp.rds")

ph <- ggplot() +
  annotate(geom='text', x=2, y=2, label = "RNA-FISH") +
  theme_void() +
  theme(panel.border = element_rect(color='black', fill=NA, size=5))

p <- y_genes + late_sperm + ph + plot_annotation(tag_levels = 'A', theme=theme(plot.caption = element_text(hjust=0)))


ggsave(snakemake@output[[1]], p, width = 10, height = 5, scale = 2)



