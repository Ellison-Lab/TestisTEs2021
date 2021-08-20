library(tidyverse)
library(GenomicRanges)
library(arrow)
library(jsonlite)

# -----------------------------------------------------------
# reviewer 1 count request
# -----------------------------------------------------------

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

lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")


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

# -----------------------------------------------------------
# pigment and epithelial cell tes
# -----------------------------------------------------------


te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv')

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr", format='arrow')

df <- map_df(w1118.obs %>% collect() %>% pull(clusters) %>% unique() %>% as.list %>% set_names(.,.),
             ~{filter(w1118.expr, clusters == . & gene_id %in% unique(te.lookup$merged_te)) %>% collect()}) %>%
  dplyr::select(index, gene_id, expression) %>%
  mutate(expression = exp(expression) - 1) %>%
  left_join(collect(w1118.obs), by=c(index='X1')) %>%
  left_join(rename.table) %>%
  group_by(clusters.rename,gene_id) %>%
  filter(sum(expression > 0) >= n()/2) %>% # TEs expressed in at least a third of the cluster
  summarise(expression = sum(expression)) %>%
  ungroup()

df %>% 
  filter(str_detect(clusters.rename,"Pig") | str_detect(clusters.rename,"Termin")) %>%
  mutate(clusters.rename=droplevels(clusters.rename)) %>% 
  split(.,.$clusters.rename) %>% 
  map(~pull(.,gene_id)) %>%
  {.[["9/Pigment"]][.[["9/Pigment"]] %in% .[["10/TerminalEpithelial"]]]}

df %>% 
  filter(str_detect(clusters.rename,"1/Spermatogonia")) %>%
  pull(gene_id)

df %>% 
  filter(str_detect(clusters.rename,"Pigment")) %>%
  pull(gene_id)

df %>% 
  filter(str_detect(clusters.rename,"Terminal")) %>%
  pull(gene_id)

df %>% 
  filter(str_detect(clusters.rename,"Cyst")) %>%
  pull(gene_id) %>%
  unique()

df %>% 
  filter(str_detect(clusters.rename,"3/Spermatocyte")) %>%
  pull(gene_id) %>%
  unique()


# -
# Diffexpression of important spermatocyte program regulators
# -
read_tsv("results/finalized/larval-w1118-testes.diffs.tsv.gz") -> x
x %>% filter(clusters.rename == "3/Spermatocyte") %>% 
  filter(gene_id %in% c("FBgn0011569","FBgn0002842","FBgn0011206","FBgn0004372",
                        "FBgn0033770","	FBgn0032473","FBgn0001313","FBgn0267432",
                        "FBgn0267433","FBgn0051702","FBgn0051703",
                        "FBgn0001313","FBgn0046323","FBgn0046697","FBgn0267432","FBgn0267433","FBgn0267592"))
