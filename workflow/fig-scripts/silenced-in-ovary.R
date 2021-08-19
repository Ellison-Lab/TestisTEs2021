library(tidyverse)
library(arrow)
library(fgsea)
library(gt)
library(ggpubr)
library(tidyverse)
library(arrow)
library(ragg)
library(rtracklayer)
library(jsonlite)
library(broom)


source("workflow/fig-scripts/theme.R")

## Get TEP-TEs

lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt") %>%
  dplyr::select(gene_id,merged_te)

geps <- open_dataset('results/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

top_mods <- geps %>%
  filter(qval < 0.005) %>%
  group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
  arrange(-n_tes)

tep_tes <- geps %>%
  filter(qval < 0.005) %>%
  filter(module == top_mods$module[1]) %>%
  left_join(lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  pull(X1) %>%
  unique()

te_mods <- geps %>%
  filter(qval < 0.005) %>%
  #filter(module == top_mods$module[1]) %>%
  left_join(lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  dplyr::select(X1,module) %>%
  split(.,.$module) %>%
  map(~pull(.,X1)) %>%
  map(unique)

res <- read_tsv("results/finalized/pirna_kd_rnaseq/pirna_kd_vs_control.res.tsv")

res <- res %>% 
  mutate(comparison=str_remove(comparison,"condition_")) %>%
  mutate(comparison=str_remove(comparison,"_vs_control"))

res <- left_join(res,lookup) %>%
  mutate(feature.type = ifelse(!is.na(merged_te),"gene","TE")) %>%
  filter(!str_detect(gene_id,"[-_]LTR"))

filter(res,!str_detect(gene_id,"FBgn")) %>% filter(is.na(merged_te)) %>%
  nrow() %>%
  {.==0} %>%
  stopifnot()

res %>%
  #filter(comparison %in% c("piwi","aub")) %>%
  filter(!str_detect(gene_id,"FBgn")) %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  pull(gene_id) %>%
  unique() -> is_silenced

df <- lookup %>%
  mutate(GEP = ifelse(merged_te %in% tep_tes,"module 27","other")) %>%
  mutate(silenced.in.ovary = gene_id %in% is_silenced) %>%
  filter(!str_detect(gene_id,"[-_]LTR")) %>%
  distinct() %>%
  group_by(GEP,silenced.in.ovary) %>%
  tally() %>%
  group_by(GEP) %>%
  mutate(prop = n/sum(n)) %>%
  mutate(GEP = fct_rev(GEP))

g <- df %>%
  ggplot(aes(GEP,prop,fill=ifelse(silenced.in.ovary,"yes","no"))) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  ylab("") +
  scale_fill_grey(limits=c("yes","no"), name="Silenced in ovary?") +
  theme_gte21() +
  xlab("module")

fish <- df %>% dplyr::select(GEP,silenced.in.ovary, n) %>% spread(silenced.in.ovary,n) %>%
  arrange(desc(GEP)) %>%
  ungroup() %>%
  dplyr::select(GEP,silenced = `TRUE`,not.silenced=`FALSE`) %>%
  column_to_rownames('GEP') %>%
  fisher.test()

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 3, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

stats.export <- fish %>%
  broom::tidy() %>%
  mutate(script= "silenced-in-ovary.R") %>%
  mutate(desc = "enrichment of TEs silenced in the ovaries in the TEP") %>%
  mutate(func = "stats::fisher.test") %>%
  mutate(ci = map2_chr(conf.low,conf.high,~paste(round(.x,digits = 3),round(.y,digits=3),sep=" - "))) %>%
  mutate(comparison = "TEP vs other") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic=estimate, ci)

write_tsv(stats.export,snakemake@output[['stats']])


