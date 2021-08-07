library(tidyverse)
library(arrow)
library(ggpubr)
library(ggpointdensity)
library(jsonlite)
library(ragg)
library(VariantAnnotation)

source("workflow/fig-scripts/theme.R")

snps <- readVcfAsVRanges("results/finalized/wgs/w1118_male/snps.vcf") %>% as_tibble()

allele.lookup <- snps %>%
  dplyr::select(seqnames, pos=start, ref, alt, specificity) %>%
  filter(specificity == 'w1118_male') %>%
  mutate(ref = ifelse(ref == "N", alt, ref))

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
  arrange(seqnames, pos)

dna <- open_dataset('results/finalized/wgs/w1118_male/pileups/', format='arrow') %>% 
  collect() %>%
  right_join(allele.lookup, by=c('seqnames','pos')) %>%
  #filter(nucleotide == ref | nucleotide == alt) %>%
  filter(specificity == "w1118_male") %>%
  filter(sample == "w1118_male") %>%
  dplyr::select( nucleotide, ref, alt, seqnames,pos, count) %>%
  arrange(seqnames, pos)

rna <- rna %>%
  dplyr::select(subsample, seqnames, pos, nucleotide, count) %>%
  complete(subsample, nesting(seqnames, pos), nucleotide, fill = list(count=0)) %>%
  filter(nucleotide!="-") %>% # don't register deletions - only considering snps
  distinct() # stand in for bug fix - upstream I return multiple counts for each entry in allele lookup

df <- left_join(dna, rna, by=c('seqnames','pos','nucleotide'), suffix=c('.dna','.tx'))

df <- mutate(df, sex = ifelse(nucleotide == alt, 'male',ifelse(nucleotide == ref, "ref", "other")))

df <- df %>% 
  dplyr::select(seqnames, pos, sex, nucleotide, ref, alt, count.dna, subsample, count.tx)

# get no-solo ltrs, add gep annotation
df <- te.lookup %>% dplyr::select(merged_te, gene_id) %>% distinct() %>%
  right_join(df, by=c(gene_id = 'seqnames')) %>%
  complete(gene_id, subsample) %>%
  mutate(gep=ifelse(gene_id %in% tep_tes,'TEP','other')) %>%
  filter(!str_detect(gene_id,'[-_]LTR'))

# get absolute depths at each pos and pct for each allele at each pos
df <- df %>% 
  group_by(gene_id, subsample, gep, pos) %>%
  mutate(depth.dna = sum(count.dna), depth.tx = sum(count.tx)) %>% 
  ungroup() %>%
  mutate(pct.tx = count.tx/depth.tx, pct.dna=count.dna/depth.dna) %>%
  arrange(subsample, gene_id, pos)

# get ratio of pct allele coverage (rna/wgs)
# this is a score for over/underexpression of the allele
df <- df %>% mutate(ratio = pct.tx/pct.dna) %>% drop_na()


# ------ Allelic expression at most male-biased sites ----------------------------

df2 <- df %>%
  filter(sex == "male") %>%
  group_by(gene_id, pos, nucleotide) %>%
  summarise(pct.tx = mean(pct.tx), ratio=mean(ratio)) %>%
  group_by(gene_id) %>%
  slice_max(ratio, n=1) %>%
  slice_max(pct.tx, n=1, with_ties = F) %>% #group_by(gene_id) %>% add_tally() %>% filter(n > 1)
  ungroup() %>%
  dplyr::select(gene_id, pos) %>%
  left_join(df) %>%  
  group_by(subsample, gene_id, pos, sex) %>%
  slice_max(ratio, n=1, with_ties = F) %>%
  ungroup() %>%
  mutate(gep = fct_relevel(gep, c("TEP","other")))

g1 <- df2 %>% 
  filter(sex == "male") %>% 
  group_by(subsample, gene_id, gep, pos) %>%
  summarise(ratio = mean(ratio)) %>%
  ggplot(aes(gep, ratio)) +
  geom_boxplot(fill="darkgray", outlier.shape = NA) +
  #geom_jitter() +
  stat_compare_means(label.y.npc = 0.9) +
  #stat_compare_means(method.args = list(alternative = "greater")) +
  theme_gte21() +
  theme(aspect.ratio = 1) + ylab("RNA/WGS") + xlab('')
  guides(fill=F)

g2 <- ggplot(df2, aes(gep, ratio)) +
  geom_boxplot(aes(fill=sex)) +
  #geom_jitter() +
  stat_compare_means(label.y.npc = 0.7) +
  theme_gte21() +
  theme(aspect.ratio = 1) + ylab("RNA/WGS") + xlab('') +
  scale_fill_brewer(type='qual', name='Variant class') +
  facet_grid(sex~subsample)

g3 <- df2 %>% 
  filter(sex != "male") %>% 
  group_by(subsample, gene_id, gep, pos) %>%
  summarise(ratio = mean(ratio)) %>%
  ggplot(aes(gep, ratio)) +
  geom_boxplot(fill="darkgray", outlier.shape = NA) +
  #geom_jitter() +
  stat_compare_means(label.y.npc = 0.7) +
  #stat_compare_means(method.args = list(alternative = "greater")) +
  theme_gte21() +
  theme(aspect.ratio = 1) + ylab("RNA/WGS") + xlab('') +
  scale_fill_brewer(type='qual', name='GEP')

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=12, height =12, units = 'in', scaling = 2, bitsize = 16, res = 300, background = 'transparent')
print(g3)
dev.off()


saveRDS(g1,snakemake@output[['ggp']])
saveRDS(g3,snakemake@output[['ggp2']])

write_tsv(df2,snakemake@output[['dat']])

# Export stats info -----------------------------------------------------------------------------------

raw.stats <- wilcox.test(y~x,data=g1$layers[[2]]$compute_aesthetics(data=g1$data, plot=g1)) %>%
  broom::tidy()

stats.export <- raw.stats %>%
  mutate(script= "w1118_pct_y_linked_rna_vs_wgs.R") %>%
  mutate(desc = "compare RNA/WGS allele depth") %>%
  mutate(func = "stats::wilcox.test/ggpubr::stat_compare_means") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "TEP TEs vs. other TEs") %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
