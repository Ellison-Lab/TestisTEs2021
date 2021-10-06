library(tidyverse)
library(ggpubr)
library(GenomicRanges)
library(arrow)
library(jsonlite)
library(ragg)

source("workflow/fig-scripts/theme.R")

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
  pull(gene_id) %>%
  unique()

copies.male <- open_dataset('results/finalized/wgs/w1118_male/copies/', format='arrow') %>% collect()
copies.female <- open_dataset('results/finalized/wgs/w1118_female/copies/', format='arrow') %>% collect()

dat.auto <- left_join(copies.male, copies.female, by="sequence", suffix=c('.male','.female')) %>%
  filter(!sequence  %in% tes)

# left_join(copies.male, copies.female, by="sequence", suffix=c('.male','.female')) %>%
#   filter(!sequence  %in% tes) %>%
#   dplyr::select(sequence,est.copies.male,est.copies.female) %>%
#   gather(sample,copies,est.copies.male,est.copies.female) %>%
#   ggplot(aes(sequence,copies)) +
#   geom_col() +
#   facet_wrap(~sample) +
#   theme_classic() +
#   theme(aspect.ratio = 1)
# 
# left_join(copies.male, copies.female, by="sequence", suffix=c('.male','.female')) %>%
#   filter(!sequence  %in% tes) %>%
#   dplyr::select(sequence,median.cov.male,median.cov.female) %>%
#   gather(sample,median.cov,median.cov.male,median.cov.female) %>%
#   ggplot(aes(sequence,median.cov)) +
#   geom_col() +
#   facet_wrap(~sample) +
#   theme_classic() +
#   theme(aspect.ratio = 1)

dat <- left_join(copies.male, copies.female, by="sequence", suffix=c('.male','.female')) %>%
  filter(sequence %in% tes) %>%
  mutate(tep = ifelse(sequence %in% tep_tes,'module 27','other')) %>%
  filter(!str_detect(sequence,'[-_]LTR')) %>%
  filter(est.copies.male >= 1 | est.copies.female > 1) %>%
  mutate(m.f.ratio = est.copies.male/est.copies.female) %>%
  mutate(tep = fct_relevel(tep,c("module 27","other")))

g3 <- dat %>%
  group_by(sign(log2(m.f.ratio))) %>%
  top_n(10,abs(log2(m.f.ratio))) %>%
  arrange(-m.f.ratio) %>%
  ggplot(aes(reorder(sequence,m.f.ratio),log2(m.f.ratio))) +
  geom_col(aes(fill=tep)) +
  theme_gte21() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size=rel(1))) +
  xlab("") + ylab("log2(M/F copies)") +
  scale_fill_gte21("binary", name="module")

g2 <- dat %>% 
  ggplot(aes(est.copies.male,est.copies.female, label=sequence)) +
  geom_point(aes(color=tep)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_gte21() +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~tep) +
  scale_color_gte21("binary",name="module") +
  scale_size(name="M/F copies")

g1 <- dat %>%
  ggplot(aes(tep,est.copies.male/est.copies.female, label=sequence))+
  geom_boxplot(outlier.shape = NA, fill="darkgray") +
  #geom_jitter(width = 0.1, alpha = 0.5) +
  stat_compare_means(label.y.npc = 0.9,size=7/.pt) +
  xlab('') +
  theme_gte21() +
  guides(fill=F) +
  ylab("M/F copies")

agg_png(snakemake@output[['png1']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

agg_png(snakemake@output[['png3']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g3)
dev.off()


saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])
saveRDS(g3,snakemake@output[['ggp3']])

write_tsv(dat,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

raw.stats <- wilcox.test(y~x,data=g1$layers[[2]]$compute_aesthetics(data=g1$data, plot=g1)) %>%
  broom::tidy()

stats.export <- raw.stats %>%
  mutate(script= "w1118_y_linked_copies.R") %>%
  mutate(desc = "compare M/F estimated copies") %>%
  mutate(func = "stats::wilcox.test/ggpubr::stat_compare_means") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "TEP TEs vs. other TEs") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
