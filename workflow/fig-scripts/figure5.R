library(tidyverse)

library(patchwork)

extrafont::loadfonts()

larrac <- read_rds('results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds') + theme(aspect.ratio = NULL)

tidal <- read_rds('results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds') + theme(aspect.ratio = NULL) + guides(fill=F)

w1118_copies_box <- read_rds('results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.ggp.rds') + theme(aspect.ratio = NULL)

male_expr <- read_rds('results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.ggp.rds') + theme(aspect.ratio = NULL) + guides(fill=F)

layout <-"
AB
CD
"

p <- tidal + larrac + w1118_copies_box + male_expr  + plot_annotation(tag_levels = 'A') +
  plot_layout(design=layout)



ggsave(snakemake@output[[1]], p, width = 7, height = 7, scale = 2)
