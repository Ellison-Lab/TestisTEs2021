library(tidyverse)
library(arrow)
library(jsonlite)


## diffs
read_tsv("results/finalized/larval-w1118-testes.diffs.tsv.gz") %>%
  dplyr::select(-clusters.orig, -clusters_id.orig) %>%
  dplyr::rename(clusters=clusters.rename) %>%
  write_tsv("~/Downloads/supplement-differential-expression.tsv.gz")

## tep

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format='arrow')

gep <- w1118.gep_membership %>% collect()

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

geps %>%
  filter(!str_detect(X1,"FBgn")) %>%
  filter(module == 27) %>%
  filter(qval < 0.005) %>%
  dplyr::rename(membership.score=weight) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  dplyr::select(X1,membership.score,qval,class=Class,family=repFamily) %>%
  distinct() %>%
  group_by(X1) %>% 
  slice_head(n = 1) %>%
  write_tsv("~/Downloads/supplement-mod27-tes.tsv")

geps %>%
  filter(str_detect(X1,"FBgn")) %>%
  filter(module == 27) %>%
  filter(qval < 0.005) %>%
  dplyr::rename(membership.score=weight) %>%
  write_tsv("~/Downloads/supplement-mod27-genes.tsv")


