library(tidyverse)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

n_cells <- read_rds("results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.n_cells.ggp.rds") + theme(aspect.ratio = NULL)

base_size <- 10

corr_plot <- read_rds("results/figs/comparison_with_mahadevaraju/comparison_with_mahadevaraju.ggp.rds") + 
  theme(aspect.ratio = NULL) +
  theme(legend.box.spacing = unit(0,"pt"), legend.box.margin = margin())+
  theme(axis.title = element_text(size=rel(1), margin = margin())) +
  theme(          legend.background = element_rect(linetype = 1),
                  legend.spacing = unit(base_size * 0.5, "points"),
                  legend.key = element_rect(linetype = 0),
                  legend.key.size = unit(0.5, "lines"),
                  legend.key.height = NULL,
                  legend.key.width = NULL,
                  legend.text = element_text(size = rel(0.75)),
                  legend.text.align = NULL,
                  legend.title = element_text(size = rel(0.7),  hjust = 0.5),
                  legend.title.align = NULL,
                  legend.position = "right",
                  legend.direction = "vertical",
                  legend.justification = "center")
  

polya_all <- read_rds("results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.all.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(axis.title = element_text(size=rel(1)))

polya_tes <- read_rds("results/figs/larval_scrna_vs_bulk/larval_scrna_vs_bulk.tes.ggp.rds") + 
  theme(aspect.ratio = NULL, axis.title = element_text(margin = margin())) +
  theme(axis.title = element_text(size=rel(1)))

batch <- read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.batch.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(axis.title = element_text(size=rel(1)))

dubs <- read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.dubs.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(axis.title = element_text(size=rel(1)),axis.title.y = element_text(margin = margin(r=-100, unit = "pt")))

umis <- read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.log1p_umis.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(axis.title = element_text(size=rel(1)))

n_genes <- read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.n_genes.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(axis.title = element_text(size=rel(1)),axis.title.y = element_text(margin = margin(r=-100,unit = "pt")))

mito <- read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.percent_mito.ggp.rds") + theme(aspect.ratio = NULL) + 
  ylab("% mitochondrial UMIs")+
  theme(axis.title = element_text(size=rel(1)))

phase <- read_rds("results/figs/larval_scrna_basic_qc_stats/larval_scrna_basic_qc_stats.phase.ggp.rds") + theme(aspect.ratio = NULL)+
  theme(axis.title = element_text(size=rel(1)))


layout <-"
AABBCC
AABBCC
DDEEFF
DDEEFF
GGHHII
GGHHII
JJ####
JJ####
"

p <- polya_all + polya_tes + corr_plot + n_cells + batch + dubs + umis + phase + n_genes + mito + plot_annotation(tag_levels = 'A') +  
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5)))


ggsave(snakemake@output[[1]], p, width = 20, height = 25)

saveRDS(p,file=snakemake@output[[2]])
