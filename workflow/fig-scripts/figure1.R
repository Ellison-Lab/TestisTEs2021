library(tidyverse)

library(patchwork)

caption <- "Larval w1118 testis scRNA-seq cell type annotation. A) UMAP embedding of cells. Each cell is colored by its parent cluster, which corresponds to a cell type found in the larval testes. B) Dotplot shows average expression level (color) and the percentage of each cluster's cells expressing a given cell type marker."

caption <- str_wrap(caption)

umap <- read_rds('results/figs/intro_larval_umap/intro_larval_umap.ggp.rds') + coord_fixed()

markers <- read_rds('results/figs/larval_marker_expression/larval_marker_expression.ggp.rds')

p <- umap + markers + plot_annotation(tag_levels = 'A', caption = caption, theme=theme(plot.caption = element_text(hjust=0)))

ggsave(snakemake@output[[1]], p, width = 10, height = 5, scale = 2)


