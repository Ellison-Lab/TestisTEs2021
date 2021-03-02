library(tidyverse)
library(patchwork)
library(png)
library(grid)
library(magick)

img <- image_read('results/figs/larval_te_heatmap/larval_te_heatmap.png', density = 300)

g <- image_ggplot(img)

caption <- "TE expression in larval testis. A) Each TE expressed in larval testis is plotted on the y-axis while cells grouped by assigned cell type are plotted on the x-axis. Cells are colored by normalized, transformed, and scaled UMI counts to visualize each TE's expression level in each cell."

caption <- str_wrap(caption)

g <- g + labs(tag="A", caption = caption) + theme(plot.caption = element_text(hjust=0))

ggsave(snakemake@output[[1]], g, width = 5, height = 3, scale = 2)





