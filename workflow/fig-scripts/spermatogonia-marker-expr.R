library(tidyverse)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/scaled/", format='arrow')

markers2display <- c("AGO3", "vas", "bam", "aub", "p53","Dek","osa","e(y)3","can",'sa',"aly")

dat <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.numeric() %>% as.list %>% set_names(.,.),
              ~{filter(w1118.expr, clusters == . & gene_symbol %in% markers2display) %>% collect()}) %>%
  left_join(collect(w1118.obs), by=c(index='X1','cell_type'='cell_type')) %>%
  left_join(rename.table) %>%
  ungroup() %>%
  #group_by(cell_type) %>%
  mutate(ord = dense_rank(gene_symbol)) %>%
  mutate(gene_symbol = fct_relevel(gene_symbol, markers2display))

# ggplot(dat, aes(fct_reorder(gene_symbol, ord),expression)) +
#   geom_violin(aes(fill=clusters.rename),draw_quantiles = c(0.5),scale = 'width') +
#   #geom_boxplot(aes(fill=clusters.rename), outlier.shape = NA) +
#   facet_wrap(~clusters.rename, ncol = 1, scales = 'free_y') +
#   scale_fill_brewer(type='qual',palette = 8, name="") +
#   theme_classic() +
#   theme(plot.caption= element_text(hjust=0.5, face='italic', size=rel(1.2)), axis.title = element_text(size = rel(1.2)), axis.text.y = element_text(size=rel(1.2)), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=rel(1.2)),aspect.ratio = 0.1) +
#   xlab('') + guides(fill=F)

dat2 <- dat %>% 
  group_by(clusters.rename, gene_symbol) %>%
  summarize(pct.expressing = sum(expression > 0)/n(), mean.expression=mean(expression))

g <- ggplot(dat2, aes(gene_symbol, clusters.rename)) +
  geom_point(aes(size=pct.expressing, fill=mean.expression), shape=21) +
  #scale_fill_gte21(palette = "diverging", discrete = F, name = "Mean expression") +
  scale_size(range=c(0, rel(5)), name='Proportion expressing') +
  theme_gte21() +
  scale_fill_distiller(palette=7, direction=1, name="Mean expression") +
  theme(axis.text.x = element_text(angle=90, hjust=1), axis.text.y = element_text(face="italic")) +
  coord_flip() +
  xlab('') + ylab('')

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(dat2,snakemake@output[['dat']])
