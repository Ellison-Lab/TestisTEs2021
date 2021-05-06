library(tidyverse)
library(ragg)
library(ggtext)
source("workflow/fig-scripts/theme.R")


lookup <- read_tsv("results/figs/celltype_rename_table.tsv") %>%
  mutate(clusters = str_extract(clusters2,"\\d+(?=\\/.+)")) %>%
  dplyr::select(clusters, clusters.rename) %>%
  deframe()

df <- read_csv("results/finalized/x-dataset-comparison/mod_scores.csv.gz", col_types = c("ccdddc")) %>%
  mutate(clusters = ifelse(dataset=="larval",lookup[clusters],clusters))

expression <- read_csv("results/finalized/x-dataset-comparison/te_expression.csv.gz", col_types = c("ccdc"))

top_corr <- df %>%
  dplyr::select(X1, clusters, dataset) %>%
  left_join(expression,.) %>%
  group_by(feature,dataset, clusters) %>%
  summarize(mean.expr = mean(expression),.groups = "drop") %>%
  left_join(dplyr::select(filter(.,dataset=="larval" & clusters == "3/Spermatocyte"), feature,ref = mean.expr), ., by="feature") %>%
  group_by(dataset, clusters) %>%
  do(tibble(corr=cor(.$ref,.$mean.expr, method = "spearman"))) %>%
  filter(dataset =="wt")

expr_corr_df <- df %>%
  dplyr::select(X1, clusters, dataset) %>%
  left_join(expression,.) %>%
  group_by(feature,dataset, clusters) %>%
  summarize(mean.expr = mean(expression),.groups = "drop") %>%
  left_join(dplyr::select(filter(.,dataset=="larval" & clusters == "3/Spermatocyte"), feature,ref = mean.expr), ., by="feature") %>%
  left_join(top_corr,.) %>%
  #mutate(dataset = paste("Witt et al.",str_to_upper(dataset)))
  filter(dataset=="wt") %>%
  mutate(is_top_hit = corr > 0.2) %>%
  mutate(clusters = fct_reorder(clusters,corr))
  

g2 <- ggplot(expr_corr_df, aes(ref, mean.expr,color=is_top_hit)) +
  geom_point(size=1) + 
  facet_wrap(~reorder(clusters,-corr), scales="free", strip.position = "left") +
  ggpubr::stat_cor(color="black",method = "spearman",size=rel(1.5)) +
  guides(color=F) +
  theme_gte21() + 
  theme(aspect.ratio = NULL, strip.placement = "outside") +
  scale_color_gte21("binary",reverse = T) +
  xlab("Expression: L3 *w1118* 3/Spermatocyte") +
  ylab("") +
  geom_smooth(method = "lm",se = F) +
  ggtitle("Comparison with Witt et al. Wild Strain") +
  theme(axis.title.x = element_markdown())
  

agg_png(snakemake@output[['png_tes']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

saveRDS(g2,snakemake@output[['ggp_tes']])

write_tsv(expr_corr_df,snakemake@output[['dat_tes']])  
