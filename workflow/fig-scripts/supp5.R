library(tidyverse)
library(png)
library(grid)
library(magick)

extrafont::loadfonts(quiet=TRUE)

library(patchwork)

source("workflow/fig-scripts/theme.R")

tpaf <- read_rds("results/figs/larval_y_gene_dotplot/larval_y_gene_dotplot.tpaf.ggp.rds") + theme(aspect.ratio = NULL) +
  theme(legend.position = "right", legend.box = "vertical", legend.direction = "vertical") +
  guides(fill=guide_legend(title="log-norm UMIs")) +
  guides(size=guide_legend(title="Prop. expressing")) +
  theme(axis.text.x = element_text(angle = 25, hjust=1, vjust=1))
  
# late_sperm <- read_rds("results/figs/larval_later_sperm_marker_umis/larval_later_sperm_marker_umis.ggp.rds") + 
#   theme(aspect.ratio = NULL) + facet_wrap(~gene_symbol, ncol = 1) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

enr <- read_rds("results/figs/larval_tep_go_enrichment/larval_tep_go_enrichment.ggp.rds") + 
  theme(axis.title.x = element_text(size=rel(0.5)), aspect.ratio = NULL) +
  #scale_y_discrete(position = "right", label = function(x) str_trunc(x, 25)) +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  geom_text(aes(x=0.25, label=Term), size=rel(2), hjust=0)


tep_y_at_least_1 <- read_rds("results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.at_least_1.ggp.rds") + 
  theme(axis.title = element_text(size=rel(1)), axis.title.y = element_text(size=rel(1)))

y_gene_enr_datalist <- read_rds("results/figs/larval_y_gene_enr/larval_y_gene_enr.ggp.rds")

pval <- y_gene_enr_datalist[[2]] %>% round(.,6) %>% paste0("P=",.) %>% str_wrap(.,width = 20)

y_gene_enr <- y_gene_enr_datalist[[1]] + theme(aspect.ratio = NULL) +
  #geom_text(aes(x=2, y=0.75, label=pval), size=rel(1.5)) +
  theme(axis.title.x = element_text(size=rel(0.5))) +
  theme(legend.title = element_blank(), legend.text = element_text(size=rel(0.5)))

#y_copies <- read_rds("results/figs/w1118_y_linked_copies/w1118_y_linked_copies.2.ggp.rds") + 
#  theme(aspect.ratio = NULL) +
#  theme(axis.title = element_text(size=rel(1))) +
#  theme(legend.key.size = unit(0.1, "pt"), legend.text = element_text(size=rel(0.5)),
#        legend.title = element_text(size=rel(0.5)), 
#        legend.box.spacing = unit(1, "pt"), strip.text = element_text(size=rel(1))) + guides(color=F, size=F) + 
#  theme(aspect.ratio = 1) +
#  facet_wrap(~tep, ncol = 2) +
#  xlab("Male copies") + ylab("Female copies") +
#  theme(axis.text = element_text(size=rel(0.7)))

y_copies_2 <- read_rds("results/figs/w1118_y_linked_copies/w1118_y_linked_copies.3.ggp.rds") + theme(aspect.ratio = NULL) +
  theme(legend.position = "right") +
  #theme(legend.position = c(0.94,0.21), legend.title = element_blank(), legend.margin = margin(1,1,1,1), legend.text = element_text(size=rel(1))) +
  theme(axis.title.y = element_text(size=rel(0.5), margin = margin()), axis.text.y = element_text(margin = margin()), axis.text.x = element_text(face="italic"))

y_linked_expr <- read_rds("results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.box.paired.ggp.rds") + 
  theme(aspect.ratio = NULL) +
  theme(axis.title.y = element_text(size=rel(0.5)))

#chr4 <- read_rds("results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.4.ggp.rds") +
#  theme(axis.title.y = element_text(size=rel(0.5))) +
#  theme(aspect.ratio = NULL)


layout <-"
AAAAA
AAAAA
AAAAA
BBBCC
BBBCC
BBBDD
BBBDD
EE###
EE###
"


p <- tpaf +  enr + y_gene_enr + tep_y_at_least_1 + #y_copies_2 + 
  y_linked_expr + #chr4 +
  plot_annotation(tag_levels = 'A') + 
  plot_layout(design=layout) &
  theme(plot.tag = element_text(face = 'bold', size=rel(1.5))) &
  theme(text=element_text(size=unit(7,"pt")))

ggsave(snakemake@output[[1]], p, width = 10, height = 12)

saveRDS(p,file=snakemake@output[[2]])
