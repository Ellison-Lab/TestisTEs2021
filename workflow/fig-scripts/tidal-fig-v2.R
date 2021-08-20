library(tidyverse)
library(jsonlite)
library(GenomicRanges)
library(arrow)
library(tidyverse)
library(arrow)
library(ragg)
library(ggpubr)


source("workflow/fig-scripts/theme.R")
options(stringsAsFactors = F)



import_tidal <- function(x) {
  z <- paste0(x,"/*result/*Inserts_Annotated.txt")
  map_df(Sys.glob(z), read_tsv) %>%
    dplyr::select(Chr,TE, loci_code) %>%
    distinct()
}

ins <- import_tidal("subworkflows/gte21-tidal-xa/results/tidal/data/DGRP/")

te_lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")%>%
  dplyr::select(merged_te, te=gene_id) %>%
  distinct()

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

w1118.gep_membership <- open_dataset("results/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

tep.name <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  left_join(te_lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  distinct() %>%
  group_by(module) %>%
  summarize(n_tes=n()) %>%
  arrange(-n_tes) %>%
  head(1) %>%
  pull(module) %>% as.character()

tep.members <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  filter(as.character(module)==tep.name) %>%
  pull(X1)

lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt") %>% #read_tsv(snakemake@input[['lookup']]) %>%
  #dplyr::select(merged_te,Flybase_name,Blast_candidate) %>%
  mutate(Flybase_name = ifelse(is.na(Flybase_name),Blast_candidate,Flybase_name)) %>%
  dplyr::select(merged_te,Flybase_name) %>%
  distinct()

# -------
# Mappable size
# 
# https://deeptools.readthedocs.io/en/develop/content/help_faq.html#use-gem
#
# > csplit -z gem_mappability_dm6_100mer.mappability "/~chr/" {*}
# > grep "chr4" xx*
# xx05:~chr4

# get sizes of chroms
sizes <- tibble(Chr = c("chrX","chr4","chr2L","chr2R","chr3L","chr3R"),Chr_len = c(20889091,1129589,21514562,20971715,23898827,28181034))
  

# sum autosomes, no chrom4 because ancient sex chrom
sizes <- sizes %>%
  #filter(Chr!='chr4') %>%
  mutate(Chr_type = ifelse(!Chr %in% c("chrX",'chrY','chr4'),'autosome',Chr)) %>%
  group_by(Chr_type) %>%
  summarize(Chr_len = sum(Chr_len))

# TART/A/B/C appear to be inconsistently reconciled by TIDAL annotation.
# I will make sure that my top-te mod remains correct by manually changing
# TART-A's REpbase name to TART-A in the lookup.

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TART-A','TART-A',Flybase_name))

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TART-B','TART-B',Flybase_name))

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TART-C','TART-C',Flybase_name))

#rename of stalker3t to stalker3: https://github.com/bergmanlab/transposons/blob/master/current/transposon_sequence_set.readme.txt

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'Stalker3','Stalker3T',Flybase_name))

lookup <- lookup %>%
  mutate(Flybase_name = ifelse(merged_te == 'TLD2','TLD2_LTR',Flybase_name))

ins <- ins %>%
  mutate(TE = ifelse(TE=="jockey","Jockey",TE))

# only include discoverable TEs (ie TEs that remain in our scRNA dataset)
ins <- ins %>% filter(TE %in% lookup$Flybase_name)


# only consider DGRP
# do not consider 4, which may have sex chrom origins
#sex_auto_ratio_df <- ins %>%
#  filter(Chr!='chr4')

# count insertions
sex_auto_ratio_df <- group_by(ins, TE, Chr) %>%
  tally() %>%
  ungroup() %>%
  complete(Chr,TE,fill=list(n=0)) %>%
  mutate_at(vars(Chr),as.character)

# normalize to length
sex_auto_ratio_df <- mutate(sex_auto_ratio_df, Chr_type = ifelse(!Chr %in% c("chrX",'chrY',"chr4"),'autosome',Chr)) %>%
  group_by(TE,Chr_type) %>%
  summarize(n=sum(n)) %>%
  left_join(sizes) %>%
  mutate(ins_per_mb = n/(Chr_len/1e6)) %>%
  ungroup() %>%
  dplyr::select(-Chr_len,-n) %>%
  spread(Chr_type,ins_per_mb)

sex_auto_ratio_df2 <- sex_auto_ratio_df %>%
  #filter(chrY > 0) %>%
  filter((autosome) > 0 & (chrX) > 0) %>%
  mutate(chrX_Auto_ratio = (chrX)/(autosome)) %>%
  mutate(chr4_Auto_ratio = (chr4)/(autosome)) %>%
  dplyr::select(TE,chrX_Auto_ratio,chr4_Auto_ratio) %>%
  mutate(is_top_mod = ifelse(TE %in% pull(filter(lookup, merged_te %in% tep.members),Flybase_name),"module 27","other")) %>%
  gather(comparison,ratio,-TE, -is_top_mod)

# ========
# plotting
# =======

g1 <- sex_auto_ratio_df2 %>%
  filter(comparison == "chrX_Auto_ratio") %>%
  ggplot(aes(is_top_mod,ratio), fill='white') +
  stat_boxplot(outlier.shape=NA, fill="darkgray") +
  ggpubr::stat_compare_means(method = 'wilcox.test',paired = F, size=rel(4),label.y = 2.5, label.x=1) +
  theme_gte21() + 
  theme(aspect.ratio = 1) +
  ylab('X/A ( per mappable MB)') +
  #scale_fill_brewer(type='qual', name='GEP') +
  coord_cartesian(ylim=c(0,3)) +
  theme(axis.title.x = element_blank())

g2 <- sex_auto_ratio_df2 %>%
  filter(comparison == "chr4_Auto_ratio") %>%
  ggplot(aes(is_top_mod,ratio), fill='white') +
  stat_boxplot(outlier.shape=NA, fill="darkgray") +
  ggpubr::stat_compare_means(method = 'wilcox.test',paired = F, size=rel(4),label.y = 2.5, label.x=1) +
  theme_gte21() + 
  theme(aspect.ratio = 1) +
  ylab('4/A (per mappable MB)') +
  #scale_fill_brewer(type='qual', name='GEP') +
  coord_cartesian(ylim=c(0,3)) +
  theme(axis.title.x = element_blank())

pct.reduction <- sex_auto_ratio_df2 %>% filter(comparison == "chrX_Auto_ratio") %>%
  group_by(is_top_mod) %>%
  summarise(ratio=median(ratio))

# ========
# export
# =======

agg_png(snakemake@output[['png1']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()


agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])

write_tsv(sex_auto_ratio_df2,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

raw.stats <- wilcox.test(y~x,data=g1$layers[[2]]$compute_aesthetics(data=g1$data, plot=g1)) %>%
  broom::tidy()

stats.export <- raw.stats %>%
  mutate(script= "tidal-fig-v2.R") %>%
  mutate(desc = "compare chrX/Autosome insertion ratio") %>%
  mutate(func = "stats::wilcox.test/ggpubr::stat_compare_means") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "TEP TEs vs. other TEs") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
