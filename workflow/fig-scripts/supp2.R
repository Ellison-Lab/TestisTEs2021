library(tidyverse)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

chimerics <- read_rds("results/figs/chimerics/chimerics.ggp.rds") + 
  theme(axis.text.x = element_text(vjust=1)) + theme(axis.text.x = element_text(face="italic"),
                                                     axis.text.y=element_text(size=unit(5,"pt")))

chimerics_supporting_reads <- read_rds("results/figs/chimerics/chimerics.supporting_reads.ggp.rds") + theme(axis.text.x = element_text(face="italic"))

full_length_box <- read_rds("results/figs/w1118_full_length_tx/w1118_full_length_tx.heat.ggp.rds")

full_length_heat <- read_rds("results/figs/w1118_full_length_tx/w1118_full_length_tx.line.ggp.rds") + theme(axis.title.x = element_blank())

raw <- read_rds("results/figs/larval_raw_te_umis/larval_raw_te_umis.ggp.rds") + theme(aspect.ratio = NULL)

normal <- read_rds("results/figs/larval_normalized_te_umis/larval_normalized_te_umis.ggp.rds") + theme(aspect.ratio = NULL)

polya_all <- read_rds("results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.all.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(text=element_text(size=unit(7,"pt")))
  #theme(axis.title = element_text(size=rel(1)))

polya_tes <- read_rds("results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.tes.ggp.rds") + 
  theme(aspect.ratio = NULL, axis.title = element_text(margin = margin())) +
  theme(text=element_text(size=unit(7,"pt")))
  #theme(axis.title = element_text(size=rel(1)))


layout <-"
AABB
AABB
CCDD
CCDD
CCDD
CCDD
EEFF
EEFF
EEGG
EEGG
EEHH
EEHH
"

p <- raw + normal + 
  polya_all + polya_tes + 
  full_length_heat + full_length_box + chimerics + chimerics_supporting_reads +  plot_annotation(tag_levels = 'A') +  
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5))) &
  theme(text=element_text(size=unit(7,"pt")))

#ggsave(snakemake@output[[1]], p, width = 20, height = 20)

saveRDS(p,file=snakemake@output[[2]])

# ------------- pg ------------------------------
library(plotgardener)
#fl <- "~/Desktop/supp1.pdf"
fl <- snakemake@output[[1]]
cairo_pdf(fl,width = 7.1, height = 10)

pageCreate(width = 180, height = 254, default.units = "mm")

plotText(label = "a", fontsize = 7, x = 0, y = 5, just = "left", default.units = "mm")
plotText(label = "b", fontsize = 7, x = 90, y = 5, just = "left", default.units = "mm")
plotText(label = "c", fontsize = 7, x = 0, y = 60, just = "left", default.units = "mm")
plotText(label = "d", fontsize = 7, x = 90, y = 60, just = "left", default.units = "mm")
plotText(label = "e", fontsize = 7, x = 0, y = 135, just = "left", default.units = "mm")
plotText(label = "f", fontsize = 7, x = 85, y = 130, just = "left", default.units = "mm")
plotText(label = "g", fontsize = 7, x = 85, y = 170, just = "left", default.units = "mm")
plotText(label = "h", fontsize = 7, x = 85, y = 207, just = "left", default.units = "mm")

plotGG(raw,x = 0,y = 7,width = 90, height = 55, default.units = "mm")
plotGG(normal,x = 85,y = 7,width = 90, height = 55, default.units = "mm")
plotGG(polya_all,x = 5,y = 50,width = 80, height = 80, default.units = "mm")
plotGG(polya_tes,x = 95,y = 50,width = 65, height = 80, default.units = "mm")
plotGG(full_length_heat,x = 5,y = 130,width = 80, height = 120, default.units = "mm")

plotGG(full_length_box,x = 85,y = 120,width = 80, height = 57, default.units = "mm")
plotGG(chimerics,x = 85,y = 170,width = 80, height = 40, default.units = "mm")
plotGG(chimerics_supporting_reads,x = 85,y = 195,width = 80, height = 65, default.units = "mm")

pageGuideHide()

dev.off()

