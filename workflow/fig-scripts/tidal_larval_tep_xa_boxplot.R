library(tidyverse)
library(arrow)
library(ragg)
library(ggpubr)


source("workflow/fig-scripts/theme.R")

df <-  open_dataset("results/finalized/larval-w1118-testes/xa_ratio/", format='arrow') %>% collect() %>% 
  filter(comparison=='chrX_Auto_ratio') %>%
  mutate(GEP = fct_relevel(GEP,c("TEP","other")))

g <- ggplot(df, aes(GEP,ratio), fill='white') +
  stat_boxplot(outlier.shape=NA, fill="darkgray") +
  ggpubr::stat_compare_means(method = 'wilcox.test',paired = F, size=rel(4),label.y = 2.5, label.x=1) +
  theme_gte21() + 
  theme(aspect.ratio = 1) +
  ylab('X/A (insertions per mb)') +
  #scale_fill_brewer(type='qual', name='GEP') +
  coord_cartesian(ylim=c(0,3)) +
  theme(axis.title.x = element_blank())

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])


# Export stats info -----------------------------------------------------------------------------------

raw.stats <- wilcox.test(y~x,data=g$layers[[2]]$compute_aesthetics(data=g$data, plot=g)) %>%
  broom::tidy()

stats.export <- raw.stats %>%
  mutate(script= "tidal_larval_te_xa_boxplot.R") %>%
  mutate(desc = "compare chrX/Autosome insertion ratio") %>%
  mutate(func = "stats::wilcox.test/ggpubr::stat_compare_means") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "TEP TEs vs. other TEs") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
