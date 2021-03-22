library(tidyverse)
library(arrow)
library(ragg)
library(rtracklayer)
library(jsonlite)

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format='arrow')
optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()
w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')
w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr/", format='arrow')

gtf <- import('~/work/TestisTpn/data/combined.fixed.gtf') %>%
  as_tibble()

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

tep.name <-  %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  dplyr::select(module, X1, repClass, repFamily) %>%
  distinct() %>%
  group_by(module) %>%
  summarize(n_tes=n()) %>%
  arrange(-n_tes) %>%
  head(1) %>%
  pull(module) %>% as.character()

geps <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect()

ylinked <- gtf %>%
  filter(seqnames == 'Y' & type == 'mRNA') %>% 
  dplyr::select(seqnames, gene_symbol, gene_id) %>%
  distinct() %>%
  filter(gene_id %in% geps$X1)

other.detected <- geps %>% filter(str_detect(X1,'FBgn')) %>% pull(X1) %>% unique() %>% length()
  
geps %>%
mutate(is.y = X1 %in% ylinked$gene_id) %>%
  group_by(module) %>%
  summarise(n.ylinked = sum(is.y), n.other  = sum(!is.y), n.ylinked.in.other = nrow(ylinked) - n.ylinked, n.other.in.other = other.detected - n.other ) %>%
  arrange(-n.ylinked) %>%
  filter(module == tep.name) %>%
  gather(x, n, n.ylinked.in.other, n.other.in.other)
# TODO
