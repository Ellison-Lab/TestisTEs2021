library(tidyverse)

library(patchwork)

umap <- read_rds('results/figs/larval_tep_usage_umap/larval_tep_usage_umap.ggp.rds') + coord_fixed()

barch <-  read_rds('results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.ggp.rds')  + theme(aspect.ratio = 0.5)

piech <- read_rds('results/figs/larval_tep_pie/larval_tep_pie.ggp.rds')


p <- umap + barch / piech + plot_annotation(tag_levels = 'A')



ggsave(snakemake@output[[1]], p, width = 6, height = 3, scale = 2)
