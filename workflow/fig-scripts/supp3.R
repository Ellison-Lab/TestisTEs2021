library(tidyverse)
library(ggtext)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

sizes <- read_rds("results/figs/larval_gep_size_histogram/larval_gep_size_histogram.ggp.rds") + 
  theme(aspect.ratio = NULL, axis.title = element_text(size=rel(1))) +
  theme(plot.margin = margin(0,10,0,10)) +
  theme(aspect.ratio = NULL, strip.text = element_text(size=rel(0.7))) +
  theme(axis.title = element_text(size=rel(1)), axis.text = element_text(size=rel(0.5)), axis.title.y = element_text(size=rel(0.5)), axis.ticks.length.x =  unit(1,"mm"))

enr <- read_rds("results/figs/larval_ica_optimized_enr/larval_ica_optimized_enr.ggp.rds") + 
  theme(aspect.ratio = NULL, axis.title = element_text(size=rel(1))) +
  ylab("% GEPs with unique enrichment") +
  theme(plot.margin = margin(0,10,0,10))

heats <- read_rds("results/figs/ica_reproducibility/ica_reproducibility.pval.ggp.rds") + 
  theme(aspect.ratio = NULL, strip.text = element_text(size=rel(0.7))) +
  theme(axis.title = element_text(size=rel(1)), axis.text = element_text(size=rel(0.5))) +
  theme(legend.title = element_text(size=rel(0.5))) +
  theme(aspect.ratio = 0.5, plot.caption = element_text(hjust=0))


#te_expr_scores <- read_rds("results/figs/all_dataset_tep_scores/all_dataset_tep_scores.tes.ggp.rds") + 
#  theme(aspect.ratio = NULL, axis.title.x = element_markdown(size=rel(0.5)),strip.text = element_text(size=rel(0.7))) +
#  theme(plot.title = element_text(size=rel(0.5)))


#layout <-"
#AAAAAA
#BBBCCC
#DDDDDD
#DDDDDD
#DDDDDD
#"

layout <-"
AAAAAA
BBBCCC
DDDDDD
DDDDDD
DDDDDD
"

p <- heats + sizes + enr + #te_expr_scores + 
  plot_annotation(tag_levels = 'a') + 
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5))) &
  theme(text=element_text(size=unit(7,"pt")))

ggsave(snakemake@output[[1]], p, width = 6, height = 6)

saveRDS(p,file=snakemake@output[[2]])
