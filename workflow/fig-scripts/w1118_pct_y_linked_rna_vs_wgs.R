library(tidyverse)
library(arrow)
library(ggpubr)
library(ggpointdensity)
library(jsonlite)
library(ragg)

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
  spread(.,sex, depth, fill = 0) %>%
  mutate(.,pct.male=male/(unknown + male), ratio.male=male/unknown)

dna <- open_dataset('results/finalized/wgs/w1118_male/snp_depths/', format='arrow') %>% 
  collect() %>%
  spread(.,sex, depth, fill = 0) %>%
  mutate(.,pct.male=male/(unknown + male), ratio.male=male/unknown)

df <- left_join(dna, rna, by=c('seqnames','pos'), suffix=c('.dna','.tx')) %>%
  group_by(subsample, sample.tx, seqnames) %>% 
  summarise(pct.male.dna=mean(pct.male.dna), pct.male.tx=mean(pct.male.tx),
            ratio.male.dna=mean(ratio.male.dna), ratio.male.tx=mean(ratio.male.tx), .groups = 'drop')

df <- te.lookup %>% dplyr::select(merged_te, gene_id) %>% distinct() %>%
  left_join(df, by=c(gene_id = 'seqnames')) %>%
  complete(gene_id, nesting(sample.tx, subsample)) %>%
  filter(!is.na(sample.tx)) %>%
  mutate(gep=ifelse(gene_id %in% tep_tes,'tep','other')) %>%
  mutate_at(vars(c('pct.male.dna','pct.male.tx','ratio.male.dna','ratio.male.tx')), replace_na, 0) %>%
  filter(!str_detect(gene_id,'[-_]LTR'))


g0 <- df %>%
  ggplot(aes(gep, pct.male.tx/pct.male.dna)) +
  geom_boxplot(aes(fill=gep)) +
  stat_compare_means() +
  theme_classic() +
  theme(aspect.ratio = 1) + ylab("% Male Depth (Male Allele/DNA)") + xlab('') +
  scale_fill_brewer(type='qual', name='GEP') +
  facet_wrap(~subsample)


g <- ggplot(df, aes(pct.male.dna,pct.male.tx, label=gene_id)) +
  geom_point(aes(color=gep)) +
  facet_wrap(~gep, scales = 'free') +
  stat_cor(method = 'spearman') + 
  theme_classic() +
  #ggrepel::geom_text_repel(force = 2) +
  geom_abline(slope=1, intercept=0, color='lightgray', linetype='dashed') +
  theme(aspect.ratio = 1) +
  scale_color_brewer(type='qual', name='GEP') +
  scale_x_continuous(labels=scales::percent) +
  scale_y_continuous(labels=scales::percent) +
  xlab('WGS mean male SNP depth') +
  ylab('RNA-seq mean male SNP depth')  +
  guides(color=F) +
  facet_grid(gep~subsample)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

agg_png(snakemake@output[['png0']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g0)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
saveRDS(g0,snakemake@output[['ggp0']])
write_tsv(df,snakemake@output[['dat']])
