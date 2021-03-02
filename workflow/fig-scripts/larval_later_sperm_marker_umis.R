library(tidyverse)
library(arrow)
library(ragg)

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr", format='arrow')

post.mei <- c('p-cup','wa-cup', 'twe','Dpy-30L2')

late.sperm <- tibble(gene_symbol = post.mei, group = 'late spermatocyte')

tmac <- c('aly','wuc','tomb') %>% tibble(gene_symbol = ., group='tMAC')
ttaf <- c('sa') %>% tibble(gene_symbol = ., group='tTAF')
tbrd <- c('tbrd-1','tbrd-2') %>% tibble(gene_symbol = ., group='tBRD')
tplus <- c('tplus3a','tplus3b') %>% tibble(gene_symbol = ., group='tPAF')

male.meiosis1.associated <- bind_rows(tmac, ttaf, tbrd, tplus)

markers <- bind_rows(male.meiosis1.associated, late.sperm)

df <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.list %>% set_names(.,.),
             ~{filter(w1118.expr, clusters == . & gene_symbol %in% markers$gene_symbol) %>% collect()}) %>%
  dplyr::select(index, gene_symbol, expression) %>%
  mutate(expression = exp(expression) - 1) %>%
  left_join(collect(w1118.obs), by=c(index='X1')) %>%
  distinct() %>%
  filter(str_detect(clusters2, 'Sperm'))

df <- left_join(markers, df)

df <- df %>%
  left_join(rename.table)  

g <- ggplot(df, aes(clusters.rename,expression)) +
  geom_violin(aes(fill=clusters.rename),draw_quantiles = c(0.5),scale = 'width') +
  theme_classic() +
  xlab("") + ylab('log2(Normalized UMIs + 1)') +
  guides(fill=F) +
  scale_fill_brewer(type='qual', palette = 8) +
  theme(aspect.ratio = 0.3) +
  theme(plot.caption= element_text(hjust=0.5, face='italic', size=rel(1.2)),
        axis.title = element_text(size = rel(1.2)), 
        axis.text.y = element_text(size=rel(1)),
        axis.text.x = element_text(size=rel(1.5), angle=90, hjust=1, vjust=0.5)) +
  facet_wrap(~gene_symbol, ncol = 2, scales='free_y')

agg_png(snakemake@output[['png']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
