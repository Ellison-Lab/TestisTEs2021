library(tidyverse)
library(ggtext)
library(patchwork)

extrafont::loadfonts()

larrac <- read_rds('results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds') + theme(aspect.ratio = NULL) +
  theme(axis.title.y = element_text(size=rel(0.5)), axis.title.x=element_blank())

tidal <- read_rds('results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds') + theme(aspect.ratio = NULL) + guides(fill=F) +
  theme(axis.title.y = element_text(size=rel(0.5)))

w1118_copies_box <- read_rds('results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.ggp.rds') + theme(aspect.ratio = NULL) +
  theme(axis.title.y = element_text(size=rel(0.5)))

male_expr <- read_rds('results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.ggp.rds') + theme(aspect.ratio = NULL) + guides(fill=F) +
  theme(axis.title.y = element_text(size=rel(0.5)))

pirna <- read_rds("results/figs/larval_pirna_expression/larval_pirna_expression.ggp.rds") +  theme(aspect.ratio = NULL) +
  theme(axis.text.y = element_markdown())

layout <-"
ABE
CDE
"

p <- tidal + larrac + w1118_copies_box + male_expr  +  pirna + plot_annotation(tag_levels = 'A') +
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))

ggsave(snakemake@output[[1]], p, width = 7, height = 4, scale = 2)

saveRDS(p,file=snakemake@output[[2]])
