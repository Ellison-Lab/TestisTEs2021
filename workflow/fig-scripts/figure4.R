library(tidyverse)
library(ggtext)
library(patchwork)
library(magick)

extrafont::loadfonts()

y_genes <- read_rds('results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.ggp.rds') +
  theme(aspect.ratio = NULL, strip.text = element_blank())

tes <- read_rds("results/figs/larval_fish_candidate_umis/larval_fish_candidate_umis.ggp.rds") +
  theme(axis.title.y=element_text(size=rel(0.5)))

ph <- image_read("/media/mlawlor/T7/microscopy_figs/210215_accord2_calfluor610_eachm_quasar670_3p4-2.slices_1_1.representative.png") %>% image_ggplot()


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

