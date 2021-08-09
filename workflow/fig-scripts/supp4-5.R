library(tidyverse)
library(png)
library(grid)
library(magick)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

flam<- read_rds("results/figs/flam_tep_enrichment/flam_tep_enrichment.ggp.rds") + ylab("") + theme(plot.margin = margin()) + xlab("")

ovary <- read_rds("results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.ggp.rds") + theme(plot.margin = margin()) + xlab("")

layout <-"
AABB
AABB
"

p <- flam + ovary + plot_annotation(tag_levels = 'A') +  
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))

ggsave(snakemake@output[[1]], p, width = 8, height = 4)

saveRDS(p,file=snakemake@output[[2]])
