library(tidyverse)
library(rtracklayer)
library(jsonlite)
library(arrow)
library(ragg)

gtf <- import('~/work/TestisTpn/data/combined.fixed.gtf') %>%
  as_tibble()

te.lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

geps <- open_dataset('results/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

top_mods <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
  arrange(-n_tes)

tep_tes <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  filter(module == top_mods$module[1]) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  pull(gene_id) %>%
  unique()

tep.name <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  dplyr::select(module, X1, repClass, repFamily) %>%
  distinct() %>%
  group_by(module) %>%
  summarize(n_tes=n()) %>%
  arrange(-n_tes) %>%
  head(1) %>%
  pull(module) %>% as.character()

feature_types <- c('mRNA','ncRNA','snRNA','pre_miRNA','miRNA','pseudogene','snoRNA')
feature_df <- gtf %>% 
  mutate(type=as.character(type)) %>%
  filter(type %in% feature_types) %>%
  dplyr::select(gene_id,type) %>%
  distinct()

df <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  left_join(feature_df, by=c('X1'='gene_id')) %>%
  filter(module==tep.name) %>%
  mutate(type=ifelse(!str_detect(X1,'FBgn'),'TE',type)) %>%
  count(type) %>%
  mutate(pct=n/sum(n))

g <- ggplot(df, aes("",pct,fill=type)) +
  geom_bar(color='white',stat='identity') +
  coord_polar('y',start=0) +
  scale_fill_brewer(type='qual',palette=6,direction = -1, name='') +
  theme_void() +
  theme(aspect.ratio = 1, legend.position='bottom') +
  xlab('') + ylab('') +
  geom_text(aes(label = scales::percent(round(pct,3))), position = position_stack(vjust = 0.5))


agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])
