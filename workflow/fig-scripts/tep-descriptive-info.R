library(tidyverse)
library(GenomicRanges)
library(arrow)
library(jsonlite)

te.lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

tes <- te.lookup$gene_id %>% unique

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

geps <- open_dataset('results/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

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

top_mods <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
  arrange(-n_tes)

tep_tes <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  filter(module == tep.name) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  pull(X1) %>%
  unique()


# -------------
#
# ------------
lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")


read_tsv("results/figs/flam_tep_enrichment/flam_tep_enrichment.dat.tsv")

read_tsv("results/figs/ovary_silenced_enrichment/ovary_silenced_enrichment.dat.tsv")

# -----------------------------------------------------------
# reviewer 1 count request
# -----------------------------------------------------------

# via wgs
wgs_y_tep_tes <- read_tsv("results/figs/w1118_y_linked_copies/w1118_y_linked_copies.dat.tsv") %>%
  filter(m.f.ratio > 1) %>%
  left_join(lookup, by=c(sequence = 'gene_id')) %>%
  filter(tep == "TEP")


# via larrac
larrac_y_tep_tes <- read_tsv("results/figs/larracuente_y_ins_barchart/larracuente_y_ins_barchart.dat.tsv") %>%
  filter(str_detect(chrom,"Y")) %>%
  filter(prop.missing < 0.1) %>%
  dplyr::select(te,GEP) %>%
  arrange(te) %>%
  distinct() %>%
  left_join(lookup, by=c(te='gene_id')) %>%
  filter(GEP=="TEP")

# via flam
# cat dmel-all-chromosome-r6.22.fasta.out | rmsk2bed > dmel-all-chromosome-r6.22.fasta.out.bed
rpts <- read_tsv("resources/dmel-all-chromosome-r6.22.fasta.out.bed",col_names = c("chr","start","end","te","score","strand"),col_types = "cnncnc")

flam <- "X:21631891-21790731" %>% GRanges()

flam_tep_tes <- rpts %>% filter(chr == "X") %>% GRanges() %>%
  subsetByOverlaps(flam) %>%
  .$te %>%
  unique() %>%
  str_remove("-int") %>%
  tibble(te=.) %>% left_join(lookup, by=c("te"="gene_id")) %>%
  filter(merged_te %in% tep_tes) %>%
  dplyr::select(merged_te) %>%
  distinct()


unique(as.vector(unlist(c(flam_tep_tes, pull(larrac_y_tep_tes,merged_te), pull(wgs_y_tep_tes, merged_te)))))
