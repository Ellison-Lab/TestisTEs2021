library(tidyverse)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

chimerics <- read_rds("results/figs/chimerics/chimerics.ggp.rds") + theme(axis.text.x = element_text(vjust=1)) + theme(axis.text.x = element_text(face="italic"))

chimerics_supporting_reads <- read_rds("results/figs/chimerics/chimerics.supporting_reads.ggp.rds") + theme(axis.text.x = element_text(face="italic"))

full_length_box <- read_rds("results/figs/w1118_full_length_tx/w1118_full_length_tx.heat.ggp.rds")

full_length_heat <- read_rds("results/figs/w1118_full_length_tx/w1118_full_length_tx.line.ggp.rds") + theme(axis.title.x = element_blank())

raw <- read_rds("results/figs/larval_raw_te_umis/larval_raw_te_umis.ggp.rds") + theme(aspect.ratio = NULL)

normal <- read_rds("results/figs/larval_normalized_te_umis/larval_normalized_te_umis.ggp.rds") + theme(aspect.ratio = NULL)

layout <-"
AABB
AABB
CCDD
CCDD
CCEE
CCEE
CCFF
CCFF
"

p <- raw + normal + full_length_heat + full_length_box + chimerics + chimerics_supporting_reads +  plot_annotation(tag_levels = 'A') +  
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))

ggsave(snakemake@output[[1]], p, width = 20, height = 20)

saveRDS(p,file=snakemake@output[[2]])
