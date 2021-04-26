library(tidyverse)
library(patchwork)
library(png)
library(grid)
library(magick)

extrafont::loadfonts()

img <- image_read('results/figs/larval_te_heatmap/larval_te_heatmap.png', density = 300)

g <- image_ggplot(img)

g <- g + labs(tag="A")

ggsave(snakemake@output[[1]], g, width = 5, height = 3, scale = 2)





