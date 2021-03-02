library(tidyverse)
library(arrow)
library(ragg)


w1118.grid <- open_dataset("results/finalized/larval-w1118-testes/grid_enr/", format='arrow')

df <- w1118.grid %>%
  collect() %>%
  mutate_at(vars('cov','pct_unique_term','pct_enr'), scales::rescale) %>%
  mutate(score = cov * pct_unique_term)

g <- df %>%
  group_by(qval,comps) %>%
  summarize(score=mean(score)) %>%
ggplot(aes(as.factor(comps),as.factor(qval),fill=score)) +
  geom_raster(interpolate = F) +
  theme_minimal() +
  xlab("k") + ylab("qval") +
  scale_fill_viridis_c() +
  theme(aspect.ratio = 1, plot.caption = element_text(hjust=0)) +
  labs(caption = "Joint score displayed by FDR cutoff and k.") #+
  #facet_wrap(~rep)


agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
