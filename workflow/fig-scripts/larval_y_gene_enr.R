library(tidyverse)
library(arrow)
library(ragg)
library(rtracklayer)
library(jsonlite)
library(broom)

source("workflow/fig-scripts/theme.R")

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership/", format='arrow')

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

w1118.obs <- open_dataset("results/finalized/larval-w1118-testes/obs", format='arrow')

w1118.expr <- open_dataset("results/finalized/larval-w1118-testes/expr/", format='arrow')

gtf <- import('subworkflows/gte21-custom-genome/results/custom-genome/combined.fixed.gtf') %>%
  as_tibble()

rename.table <- read_tsv('results/figs/celltype_rename_table.tsv') %>%
  mutate(clusters.rename = fct_reorder(clusters.rename,as.numeric(str_extract(clusters.rename,"\\d+")))) %>%
  arrange(clusters.rename)

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

geps <- w1118.gep_membership %>%
  #
  collect()

tep.name <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  dplyr::select(module, X1, repClass, repFamily) %>%
  distinct() %>%
  group_by(module) %>%
  summarize(n_tes=n()) %>%
  arrange(-n_tes) %>%
  head(1) %>%
  pull(module) %>% as.character()

#tep_genes <- geps %>%
#  filter(qval < optimal_ica[['qval']]) %>%
#  filter(module == tep.name) %>% pull(X1)

expressed_genes <- unique(pull(collect(w1118.expr),"gene_id"))

ylinked <- gtf %>%
  filter(seqnames == 'Y' & type == 'mRNA') %>%
  dplyr::select(seqnames, gene_symbol, gene_id) %>%
  filter(gene_id %in% expressed_genes) %>%
  distinct()
  
tep <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  filter(str_detect(X1,"FBgn")) %>%
  mutate(chrom = ifelse(X1 %in% ylinked$gene_id,"Y","other")) %>%
  mutate(GEP = ifelse(module == tep.name,"TEP","non-TEP")) %>%
  filter(GEP == "TEP") %>%
  dplyr::select(X1, chrom, GEP)

non_tep <- geps %>%
  filter(str_detect(X1,"FBgn")) %>%
  mutate(chrom = ifelse(X1 %in% ylinked$gene_id,"Y","other")) %>%
  mutate(GEP = ifelse(module == tep.name,"TEP","non-TEP")) %>%
  filter(GEP == "non-TEP") %>%
  filter(!X1 %in% tep$X1) %>%
  dplyr::select(X1, chrom, GEP) %>%
  distinct()

df0 <- bind_rows(tep, non_tep)

df <- df0 %>%
  group_by(chrom, GEP) %>%
  tally() %>%
  ungroup() %>%
  mutate(GEP = ifelse(GEP == "TEP","module 27","other"))
  

g1 <- df %>%
  group_by(chrom) %>%
  mutate(pct = n/sum(n)) %>%
  ggplot(aes(chrom,pct, fill=fct_rev(GEP))) +
  geom_col(position = "stack") +
  theme_gte21() +
  scale_fill_gte21("binary",reverse = T,name="module") +
  theme(aspect.ratio = 1) +
  scale_y_continuous(labels=scales::percent) +
  ylab("") + xlab("genomic location")

fish_res <- df %>%
  spread(chrom,n) %>%
  arrange(desc(GEP)) %>%
  dplyr::select(GEP,Y, other) %>%
  column_to_rownames("GEP") %>%
  fisher.test() %>%
  tidy()

pval <-  pull(fish_res,p.value)

#pool <- geps %>%
#  filter(str_detect(X1,"FBgn")) %>%
#  pull

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 3, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

saveRDS(list(g1,pval),snakemake@output[['ggp']])
write_tsv(df0,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

stats.export <- fish_res %>%
  mutate(script= "larval_y_gene_enr.R") %>%
  mutate(desc = "Y gene enrichment in TEP") %>%
  mutate(func = "stats::fisher.test") %>%
  mutate(ci = map2_chr(conf.low,conf.high,~paste(round(.x,digits = 3),round(.y,digits=3),sep=" - "))) %>%
  mutate(comparison = "TEP vs other") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic=estimate, ci)

write_tsv(stats.export,snakemake@output[['stats']])
