library(tidyverse)

library(patchwork)

extrafont::loadfonts()

source("workflow/fig-scripts/theme.R")

umap <- read_rds('results/figs/larval_tep_usage_umap/larval_tep_usage_umap.ggp.rds') +
  theme(axis.text = element_text(size=rel(1)), axis.title = element_text(size=rel(1))) +
  theme(aspect.ratio = NULL) +
  theme(legend.position = c(0.9,0.13))


barch <-  read_rds('results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.ggp.rds')  + theme(aspect.ratio = NULL) +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=rel(0.5)), axis.title.y =  element_text(size=rel(0.5)), axis.text.y=element_text(size=rel(0.5)))

piech <- read_rds('results/figs/larval_tep_pie/larval_tep_pie.ggp.rds') +
  theme(legend.position = "bottom", legend.text = element_text(size=rel(0.7)), legend.key.size = unit(1,"line"), legend.margin = margin(t=0,unit="line")) +
  theme(plot.margin = margin()) +
  theme(legend.direction = "horizontal", legend.position = c(0.55,0.05)) +
  ggtitle("GEP-27")


layout <-"
AAAAAABBB
AAAAAACCC
AAAAAACCC
"


p <- umap + barch + piech +
  plot_annotation(tag_levels = 'A', theme=theme(plot.caption = element_text(hjust=0, family="Arial"))) +
  plot_layout(design = layout)


ggsave(snakemake@output[[1]], p,  width = 10, height = 6)
