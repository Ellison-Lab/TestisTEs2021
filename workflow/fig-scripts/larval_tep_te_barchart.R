library(tidyverse)
library(arrow)
library(ragg)
library(jsonlite)

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

te.lookup <- read_tsv('~/work/TestisTpn/data/te_id_lookup.curated.tsv.txt')

tep.name <- w1118.gep_membership %>%
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

df <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  filter(module == tep.name) %>%
  collect() %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  dplyr::select(module, X1, repClass, repFamily) %>%
  distinct() %>%
  group_by(module) %>%
  mutate(n_tes=n()) %>%
  filter(n_tes > 1) %>%
  mutate(repClass=replace_na(repClass,replace = 'other')) %>%
  mutate(module=paste0("GEP-",module))


df2 <- df %>% 
  group_by(repClass) %>%
  count(repFamily) %>%
  group_by(repClass) %>%
  mutate(pct=n/sum(n)) %>%
  mutate_if(is.character, ~replace_na(.,"other"))

g <- ggplot(df2, aes("",pct,fill=repFamily)) +
  geom_bar(color='white',stat='identity') +
  coord_polar('y',start=0) +
  scale_fill_viridis_d() +
  theme_void() +
  theme(aspect.ratio = 1, legend.position='bottom') +
  xlab('') + ylab('') +
  geom_text(aes(label = repFamily), position = position_stack(vjust = 0.5), color='orange') +
  guides(fill=F) +
  #geom_text(aes(label = scales::percent(round(pct,3))), position = position_stack(vjust = 0.5)) +
  #scale_y_continuous(expand=c(0,0)) +
    facet_wrap(~repClass)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df2,snakemake@output[['dat']])
