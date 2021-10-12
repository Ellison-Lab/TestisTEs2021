library(tidyverse)

# -------------- Fig 5 -------------------------------------

larrac <- read_rds('results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.ggp.rds') + theme(aspect.ratio = NULL) +
  theme(axis.title.y = element_text(size=rel(0.5)), axis.title.x=element_blank())

#larrac$data

tidal <- read_rds('results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.ggp.rds') + theme(aspect.ratio = NULL) + guides(fill=F) +
  theme(axis.title.y = element_text(size=rel(0.5)))

tidal$data %>% group_by(is_top_mod) %>% tally

w1118_copies_box <- read_rds('results/figs/w1118_y_linked_copies/w1118_y_linked_copies.1.ggp.rds') + theme(aspect.ratio = NULL) +
  theme(axis.title.y = element_text(size=rel(0.5)))

w1118_copies_box$data %>% group_by(tep) %>% tally

male_expr <- read_rds('results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.ggp.rds') + theme(aspect.ratio = NULL) + guides(fill=F) +
  theme(axis.title.y = element_text(size=rel(0.5)))

male_expr$data %>%
  group_by(subsample,gep) %>%
  tally()


# --------------- SUpp 11 -----------------------

tidal_4_over_a <- read_rds("results/figs/tidal_larval_tep_xa_boxplot/tidal_larval_tep_xa_boxplot.4.ggp.rds") +
  theme(text=element_text(size=unit(7,"pt")))

tidal_4_over_a$data %>% group_by(is_top_mod) %>% tally()


# ----------------supp 8 e --------------------

y_linked_expr <- read_rds("results/figs/w1118_pct_y_linked_rna_vs_wgs/w1118_pct_y_linked_rna_vs_wgs.box.paired.ggp.rds") + 
  theme(aspect.ratio = NULL) +
  theme(axis.title.y = element_text(size=rel(0.5)))

y_linked_expr$data %>%
  group_by(subsample, gep) %>% tally()
