library(tidyverse)
library(png)
library(grid)
library(magick)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

pie <- read_rds("results/figs/larval_tep_te_barchart/larval_tep_te_barchart.ggp.rds") + 
  theme(legend.position = "right", legend.direction = "vertical") +
  facet_wrap(~repClass, ncol = 1)

heat1<- image_read("results/figs/larval_gep_corr_heatmap/larval_gep_corr.membership.png", density = 300) %>%
  image_trim(fuzz=1) %>%
  image_ggplot()

heat2<- image_read("results/figs/larval_gep_corr_heatmap/larval_gep_corr.usage.png", density = 300) %>%
  image_trim(fuzz=1) %>%
  image_ggplot()


heat1 <- heat1 + theme(plot.margin = margin())

heat2 <- heat2 + theme(plot.margin = margin())

layout <-"
AABB
AABB
"

p <- heat2 + heat1  +  #pie + 
  plot_annotation(tag_levels = 'a') +  
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))

ggsave(snakemake@output[[1]], p, width = 8, height = 8)

saveRDS(p,file=snakemake@output[[2]])
