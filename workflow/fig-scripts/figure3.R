library(tidyverse)

library(patchwork)

extrafont::loadfonts()

source("workflow/fig-scripts/theme.R")

umap <- read_rds('results/figs/larval_tep_usage_umap/larval_tep_usage_umap.ggp.rds') +
  theme(axis.text = element_text(size=rel(1)), axis.title = element_text(size=rel(1))) +
  theme(aspect.ratio = NULL) +
  theme(legend.position = c(0.9,0.13)) + coord_fixed()


barch <- read_rds('results/figs/larval_all_gep_te_barchart/larval_all_gep_te_barchart.ggp.rds')  + theme(aspect.ratio = NULL) +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=rel(1)), 
        axis.title.y =  element_text(size=rel(0.5)), axis.text.y=element_text(size=rel(0.5)),
        axis.title.x = element_text(size=rel(0.5))) +
  theme(legend.position = c(0.1,0.75), legend.title = element_blank(), legend.margin = margin(1,1,1,1), legend.text = element_text(size=rel(0.5)), legend.background = element_rect(color = NA)) +
  xlab("GEP")

piech <- read_rds('results/figs/larval_tep_pie/larval_tep_pie.ggp.rds') +
  theme(legend.position = "bottom", legend.text = element_text(size=rel(0.7)), legend.key.size = unit(1,"line"), legend.margin = margin(t=0,unit="line")) +
  theme(plot.margin = margin(), plot.title = element_text(hjust=0.5, margin = margin(5,5,5,5))) +
  theme(legend.direction = "horizontal", legend.position = c(0.52,0)) +
  ggtitle("GEP-27")


layout <-"
AAAACCCCCC
BBB#CCCCCC
BBB#CCCCCC
"


p <- barch + piech + umap + 
  plot_annotation(tag_levels = 'A', theme=theme(plot.caption = element_text(hjust=0, family="Arial"))) +
  plot_layout(design = layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))


ggsave(snakemake@output[[1]], p,  width = 10, height = 6)

saveRDS(p,file=snakemake@output[[2]])
