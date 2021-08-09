library(rtracklayer)
library(tidyverse)
library(arrow)
library(tidyverse)
library(arrow)
library(ragg)
library(rtracklayer)
library(jsonlite)
library(broom)


source("workflow/fig-scripts/theme.R")

lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

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

# cat dmel-all-chromosome-r6.22.fasta.out | rmsk2bed > dmel-all-chromosome-r6.22.fasta.out.bed
rpts <- read_tsv("resources/dmel-all-chromosome-r6.22.fasta.out.bed",col_names = c("chr","start","end","te","score","strand"),col_types = "cnncnc")

flam <- "X:21631891-21790731" %>% GRanges()

flam_tes <- rpts %>% filter(chr == "X") %>% GRanges() %>%
  subsetByOverlaps(flam) %>%
  .$te %>%
  unique() %>%
  str_remove("-int") %>%
  tibble(te=.) %>% left_join(lookup, by=c("te"="gene_id")) %>%
  dplyr::select(merged_te) %>%
  distinct()

flam_tes  %>% pull(merged_te) %>%
  sort()

# construct contingency tabe
flam_in_tep <- sum(flam_tes$merged_te %in% tep_tes)
flam_not_in_tep <- sum(!flam_tes$merged_te %in% tep_tes)

n_tep <- length(tep_tes)

not_flam_in_tep <- n_tep - flam_in_tep
not_flam_not_in_tep <- length(unique(lookup$merged_te)) - flam_in_tep - flam_not_in_tep - not_flam_in_tep

cont <- matrix(c(flam_in_tep, not_flam_in_tep, flam_not_in_tep, not_flam_not_in_tep),ncol = 2,dimnames = list(c("yes","no"),c("TEP","other")))

cont_df <- cont %>% as.data.frame() %>% rownames_to_column("cluster") %>%
  gather(module, count, -cluster) %>%
  group_by(module) %>%
  mutate(prop = count/sum(count))

# plot contingency table

g <- cont_df %>%
  ggplot(aes(module,prop,fill=cluster)) +
  geom_col() +
  #scale_fill_brewer(name="Detected in flamenco locus?",palette = 6, type="qual", direction = -1) +
  scale_fill_grey(name="Detected in flamenco locus?", limits=c("yes","no")) +
  theme_gte21() +
  scale_y_continuous(labels  = scales::percent)


fish <- fisher.test(cont) %>% broom::tidy()


agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 3, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(cont_df,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

stats.export <- fish %>%
  mutate(script= "flam-enrich.R") %>%
  mutate(desc = "flam-controlled TE enrichment in TEP") %>%
  mutate(func = "stats::fisher.test") %>%
  mutate(ci = map2_chr(conf.low,conf.high,~paste(round(.x,digits = 3),round(.y,digits=3),sep=" - "))) %>%
  mutate(comparison = "TEP vs other") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic=estimate, ci)

write_tsv(stats.export,snakemake@output[['stats']])
