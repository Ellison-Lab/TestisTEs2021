library(tidyverse)
library(patchwork)
library(png)
library(grid)
library(magick)

extrafont::loadfonts()

img <- image_read('results/figs/larval_te_heatmap/larval_te_heatmap.png', density = 300) %>% 
  image_trim(fuzz = 1) #%>% 
  #image_crop(geometry = geometry_area(y_off = 10))

g <- image_ggplot(img)


g2 <- read_rds("results/figs/te_expression_by_cluster/te_expression_by_cluster.ggp.rds") + 
  theme(axis.title = element_text(size=rel(0.5))) +
  #theme(aspect.ratio = 0.5) +
  theme(axis.title.x = element_blank())

layout1 <-"
A
"

layout2 <-"
B
"


p <- g + 
  plot_annotation(tag_levels = 'A', theme=theme(plot.caption = element_text(hjust=0, family="Arial"))) +
  plot_layout(design = layout1) & #,widths = 1,heights = c(3,1)) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))


p2 <- g2 + labs(tag = "B") + theme(plot.tag = element_text(face = 'bold', family="Arial", size=rel(1.5)))

ggsave(snakemake@output[[1]], p, width = 5, height = 3, scale = 2)
ggsave(snakemake@output[[2]], p2, width = 2, height = 3, scale = 2)

saveRDS(g,snakemake@output[[3]])



