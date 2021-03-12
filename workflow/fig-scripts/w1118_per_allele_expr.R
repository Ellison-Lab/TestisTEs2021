library(tidyverse)
library(arrow)
library(ggpubr)
library(ggpointdensity)
library(jsonlite)
library(ragg)
library(VariantAnnotation)

snps <- readVcfAsVRanges("results/finalized/wgs/w1118_male/snps.vcf") %>% as_tibble()

allele.lookup <- snps %>%
  dplyr::select(seqnames, pos=start, ref, alt, specificity) %>%
  filter(specificity == 'w1118_male')

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

rna <- Sys.glob('results/finalized/w1118-testes-total-rna/rep*-depth-at-male-snps/') %>%
  set_names(.,str_extract(.,"(?<=rna\\/)rep\\d+")) %>%
  map_df(~collect(open_dataset(., format='arrow'))) %>%
  distinct() %>%
  left_join(allele.lookup, by=c('seqnames','pos')) %>%
  filter(specificity == "w1118_male") %>%
  mutate(sex = ifelse(nucleotide == alt, 'male','unknown')) %>%
  group_by(seqnames, pos, nucleotide, sex, subsample) %>% mutate(n=n()) %>% arrange(-n)
  summarise(depth = sum(`count`),.groups = 'drop') %>%
  spread(.,sex, depth, fill = 0) %>%
  mutate(.,pct.male=male/(unknown + male))

dna <- open_dataset('results/finalized/wgs/w1118_male/pileups/', format='arrow') %>% 
  collect() %>%
  left_join(allele.lookup, by=c('seqnames','pos')) %>%
  filter(specificity == "w1118_male") %>%
  mutate(sex = ifelse(nucleotide == alt, 'male','unknown')) %>%
  group_by(seqnames, pos, sex, sample) %>%
  summarise(depth = sum(`count`),.groups = 'drop') %>%
  dplyr::select(sample, seqnames,pos,sex, depth) %>%
  spread(.,sex, depth, fill = 0) %>%
  mutate(.,pct.male=male/(unknown + male)) %>%
  arrange(seqnames, pos)




