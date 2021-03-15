library(tidyverse)

library(patchwork)


larrac <- read_rds('results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds')


tidal <- read_rds('results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds')


w1118_copies_box <- read_rds('results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.ggp.rds')

male_expr <- read_rds('results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.box.paired.ggp.rds') + theme(aspect.ratio = 1)

p <- ( tidal + larrac ) / w1118_copies_box + male_expr  + plot_annotation(tag_levels = 'A')



ggsave(snakemake@output[[1]], p, width = 7, height = 4, scale = 2)
