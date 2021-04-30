library(tidyverse)
library(ggtext)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

sizes <- read_rds("results/figs/larval_gep_size_histogram/larval_gep_size_histogram.ggp.rds") + 
  theme(aspect.ratio = NULL, axis.title = element_text(size=rel(1)))


heats <- read_rds("results/figs/larval_cica_grid_heatmap/larval_cica_grid_heatmap.ggp.rds") + 
  theme(aspect.ratio = NULL, strip.text = element_text(size=rel(0.7))) +
  theme(axis.title = element_text(size=rel(1)), axis.text = element_text(size=rel(0.5))) +
  theme(legend.title = element_text(size=rel(0.5)))


te_expr_scores <- read_rds("results/figs/all_dataset_tep_scores/all_dataset_tep_scores.tes.ggp.rds") + 
  theme(aspect.ratio = NULL, axis.title.x = element_markdown(size=rel(0.5)),strip.text = element_text(size=rel(0.7))) +
  theme(plot.title = element_text(size=rel(0.7)))


layout <-"
AAAAAA
BBB###
CCCCCC
CCCCCC
CCCCCC
"

p <- heats + sizes + te_expr_scores + plot_annotation(tag_levels = 'A') + 
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))

ggsave(snakemake@output[[1]], p, width = 10, height = 12)

saveRDS(p,file=snakemake@output[[2]])
