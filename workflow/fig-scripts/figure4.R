library(tidyverse)
library(ggtext)
library(patchwork)

extrafont::loadfonts()

y_genes <- read_rds('results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds') +
  theme(aspect.ratio = NULL, strip.text = element_blank())

tes <- read_rds("results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.ggp.rds") +
  theme(axis.title.y=element_text(size=rel(0.5)))

ph <- ggplot() +
  annotate(geom='text', x=2, y=2, label = "RNA-FISH") +
  theme_void() +
  theme(panel.border = element_rect(color='black', fill=NA, size=5))


layout <-"
AAAABBB
AAAABBB
AAAABBB
AAAACCC
AAAACCC
"


p <- y_genes + tes + ph + 
  plot_annotation(tag_levels = 'A', theme=theme(plot.caption = element_text(hjust=0))) +
  plot_layout(design = layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))


ggsave(snakemake@output[[1]], p, width = 10, height = 10)

saveRDS(p,file=snakemake@output[[2]])

