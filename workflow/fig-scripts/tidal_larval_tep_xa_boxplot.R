library(tidyverse)
library(arrow)
library(ragg)

df <-  open_dataset("results/finalized/larval-w1118-testes/xa_ratio/", format='arrow') %>% collect() %>% 
  filter(comparison=='chrX_Auto_ratio')
  
g <- ggplot(df, aes(GEP,ratio), fill='white') +
  stat_boxplot(outlier.shape=NA) +
  ggpubr::stat_compare_means(aes(label = paste("p =",..p.format..)),method = 'wilcox.test',paired = F, size=rel(3)) +
  theme_classic() + 
  theme(aspect.ratio = 1) +
  ylab('X/A (insertions per mb)')

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
