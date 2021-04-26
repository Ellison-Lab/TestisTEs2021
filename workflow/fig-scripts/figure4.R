library(tidyverse)

library(patchwork)

extrafont::loadfonts()

y_genes <- read_rds('results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds') +
  theme(aspect.ratio = NULL)

eachm <- read_rds("results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.eachm.ggp.rds") + ggtitle("EAChm") + ylab("log1p(expression)") +
  theme(plot.title = element_text(face=NULL))

tes <- read_rds("results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.ggp.rds") 

ph <- ggplot() +
  annotate(geom='text', x=2, y=2, label = "RNA-FISH") +
  theme_void() +
  theme(panel.border = element_rect(color='black', fill=NA, size=5))

layout <-"
AAAABBB
AAAACCC
AAAADDD
AAAADDD
AAAADDD
"


p <- y_genes + eachm + tes + ph + 
  plot_annotation(tag_levels = 'A', theme=theme(plot.caption = element_text(hjust=0))) +
  plot_layout(design = layout)


ggsave(snakemake@output[[1]], p, width = 7, height = 7, scale = 2)



