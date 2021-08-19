library(tidyverse)
library(arrow)
library(tidyverse)
library(jsonlite)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")


# ------------------------
# import usage score data
# ------------------------

final <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format="arrow") %>% collect() %>%
  mutate(components = 90)

#final <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-ica-grid/results/gep/larval-w1118-testes/optimal/consensus-usage.csv.gz") %>%
#  read_csv() %>%
#  mutate(components = 90)

reps <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-ica-grid/results/grid-search/larval-w1118-testes/overall-*/90-components/consensus-ica.csv.gz") %>%
  set_names(.,str_extract(.,"(?=overall-)\\d+")) %>%
  map_df(read_csv,.id="rep")

# --------------------------------------------------
# get corr p vals
# --------------------------------------------------
to_corr <- reps %>% gather(module, score, -rep, -index) %>%
  left_join(final,.,by=c(X1 = "index")) %>%
  dplyr::select(rep, feature="X1", module.x, score.x="weight", module.y, score.y="score")

corr_res <- to_corr %>%
  group_by(module.x, module.y, rep) %>%
  do(broom::tidy(cor.test(.$score.x, .$score.y, method="pearson", exact=T))) %>%
  ungroup()

pv.replace <- min(corr_res$p.value[corr_res$p.value > 0])

corr_res_to_plot <- corr_res %>%
  group_by(module.x, rep) %>%
  filter(estimate > 0) %>%
  slice_max(n=1, order_by = estimate) %>%
  group_by(module.x) %>%
  mutate(max.corr = max(estimate)) %>%
  ungroup() %>%
  mutate(rank = dense_rank(-max.corr)) %>%
  mutate(p.value = ifelse(p.value == 0, pv.replace, p.value)) %>%
  mutate(rep = paste("cICA replicate:", rep))

g1 <- ggplot(corr_res_to_plot,aes(rank,-log10(p.value))) +
  #stat_summary() +
  geom_jitter(width = 0.1, size=0.25) +
  geom_smooth() +
  scale_y_log10() +
  #scale_y_continuous(breaks = seq(0,100,10)/100) +
  ylab("-log10(pval)") +
  xlab("final cICA modules") +
  theme_gte21() +
  facet_wrap(~rep)


g2 <- ggplot(corr_res_to_plot,aes(rank,estimate)) +
  #stat_summary() +
  geom_jitter(width = 0.1, size=0.25) +
  geom_smooth() +
  #scale_y_continuous(breaks = seq(0,100,10)/100) +
  ylab("Spearman's rho") +
  xlab("Modules in final set\n(ranked by Spearman's rho for best match in grid search replicates)") +
  theme_gte21() +
  facet_wrap(~rep)


corr_res_to_plot%>%
  group_by(module.x) %>%
  filter(all(estimate > 0.9)) %>%
  pull(module.x) %>% unique() %>% length()

agg_png(snakemake@output[['png1']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])

write_tsv(corr_res_to_plot,snakemake@output[['dat']])

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
# ------------------------
# import usage score data
# ------------------------

# final <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_usage/", format="arrow") %>% collect() %>%
#   mutate(components = 90)
# 
# #final <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-ica-grid/results/gep/larval-w1118-testes/optimal/consensus-usage.csv.gz") %>%
# #  read_csv() %>%
# #  mutate(components = 90)
# 
# rep1 <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-ica-grid/results/grid-search/larval-w1118-testes/overall-1/*-components/consensus-usage.csv.gz") %>%
#   set_names(.,str_extract(.,"\\d+(?=-components)")) %>%
#   map_df(read_csv,.id="components") %>%
#   mutate(X1=str_remove(X1,regex("^b"))) %>%
#   mutate(X1=str_remove_all(X1,regex("\\'")))
# 
# rep2 <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-ica-grid/results/grid-search/larval-w1118-testes/overall-2/*-components/consensus-usage.csv.gz") %>%
#   set_names(.,str_extract(.,"\\d+(?=-components)")) %>%
#   map_df(read_csv,.id="components")  %>%
#   mutate(X1=str_remove(X1,regex("^b"))) %>%
#   mutate(X1=str_remove_all(X1,regex("\\'")))
# 
# rep3 <- Sys.glob("~/amarel-scratch/TE-proj-reorg/TestisTEs2021/subworkflows/gte21-ica-grid/results/grid-search/larval-w1118-testes/overall-3/*-components/consensus-usage.csv.gz") %>%
#   set_names(.,str_extract(.,"\\d+(?=-components)")) %>%
#   map_df(read_csv,.id="components")  %>%
#   mutate(X1=str_remove(X1,regex("^b"))) %>%
#   mutate(X1=str_remove_all(X1,regex("\\'")))
# 
# 
# # --------------------------------------------------
# # get corr p vals
# # --------------------------------------------------
# to_corr <- list(rep1 = rep1, rep2=rep2, rep3=rep3) %>%
#   bind_rows(.id="rep") %>%
#   filter(components==90) %>%
#   left_join(final,.,by="X1") %>%
#   dplyr::select(rep, X1, cons.x, usage.x, usage.y, cons.y)
#   
# corr_res <- to_corr %>%
#   group_by(cons.x, cons.y, rep) %>%
#   do(broom::tidy(cor.test(.$usage.x, .$usage.y, method="spearman", exact=T))) %>%
#   ungroup()
# 
# corr_res_to_plot <- corr_res %>%
#   group_by(cons.x, rep) %>%
#   slice_max(n=1, order_by = estimate) %>%
#   group_by(cons.x) %>%
#   mutate(max.corr = max(estimate)) %>%
#   ungroup() %>%
#   arrange(-max.corr) %>%
#   mutate(rank = dense_rank(p.value))
#   
# ggplot(corr_res_to_plot,aes(rank,p.value)) +
#   #stat_summary() +
#   geom_jitter(width = 0.1, size=0.1) +
#   geom_smooth() +
#   scale_y_log10() +
#   #scale_y_continuous(breaks = seq(0,100,10)/100) +
#   ylab("|Spearman's rho| for best match from 3 grid search replicates") +
#   xlab("Modules in final set\n(ranked by best match in grid search replicates)")

