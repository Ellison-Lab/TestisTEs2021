library(tidyverse)
library(arrow)
library(ragg)



w1118.gep_usage <- open_dataset("~/finalized/larval-w1118-testes/optimal_gep_usage", format='arrow')
w1118.gep_membership <- open_dataset("~/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')


w1118.membership.cor <- w1118.gep_membership %>% 
  collect() %>%
  dplyr::select(-qval) %>%
  spread(module, weight) %>%
  column_to_rownames('X1') %>% as.matrix() %>%
  cor()

w1118.usage.cor <- w1118.gep_usage %>% 
  collect() %>%
  #group_by(X1) %>%
  #mutate(usage = scales::rescale(usage,c(-1,1))) %>%
  spread(cons, usage) %>%
  column_to_rownames('X1') %>% as.matrix() %>%
  cor()

w1118.usage.hc.geps <- w1118.usage.cor %>% {1-.} %>% as.dist() %>% hclust()


w1118.membership.hc.geps <- w1118.membership.cor %>% {1-.} %>% as.dist() %>% hclust()



agg_png(snakemake@output[['png1']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
pheatmap::pheatmap(w1118.usage.cor,
                   cluster_rows = w1118.usage.hc.geps, 
                   cluster_cols = w1118.usage.hc.geps, treeheight_row = F,treeheight_col = 10,fontsize_col = 2,
                   border_color = NA, cellwidth = 3, cellheight = 3, show_rownames = F)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
pheatmap::pheatmap(w1118.membership.cor,
                   cluster_rows = w1118.membership.hc.geps, 
                   cluster_cols = w1118.membership.hc.geps, treeheight_row = F,treeheight_col = 10,fontsize_col = 2,
                   border_color = NA, cellwidth = 3, cellheight = 3, show_rownames = F)
dev.off()
