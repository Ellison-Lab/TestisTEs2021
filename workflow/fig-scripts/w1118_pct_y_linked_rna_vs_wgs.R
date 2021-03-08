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

scrna.detected.tes <- geps %>% filter(!str_detect(X1,'FBgn')) %>% dplyr::select(X1) %>% left_join(te.lookup, by=c(X1='merged_te')) %>% pull(gene_id)

rna <- Sys.glob('results/finalized/w1118-testes-total-rna/rep*-depth-at-male-snps/') %>%
  set_names(.,str_extract(.,"(?<=rna\\/)rep\\d+")) %>%
  map_df(~collect(open_dataset(., format='arrow'))) %>%
  distinct() %>% # get rid of apparent dup rows caused by multiple male-specific snps in same pos
  spread(.,sex, depth, fill = 0) %>%
  mutate(.,pct.male=male/(unknown + male))

dna <- open_dataset('results/finalized/wgs/w1118_male/snp_depths/', format='arrow') %>% 
  collect() %>%
  spread(.,sex, depth, fill = 0) %>%
  mutate(.,pct.male=male/(unknown + male))

res2 <-inner_join(dna, rna, by=c('seqnames','pos'), suffix=c('.dna','.tx')) %>%
  group_by(sample.tx, seqnames) %>% 
  summarise(pct.male.dna=mean(pct.male.dna), pct.male.tx=mean(pct.male.tx), .groups = 'drop')

res3 <- res2%>% mutate(gep=ifelse(seqnames %in% tep_tes,'tep','other'))

df <- res3 %>%
  filter(!str_detect(seqnames,'[-_]LTR')) %>%
  filter(seqnames %in% scrna.detected.tes)

g0 <- df %>% mutate(ratio=pct.male.tx/pct.male.dna) %>% 
  ggplot(aes(gep, ratio)) +
  geom_boxplot(aes(fill=gep)) +
  stat_compare_means() +
  theme_classic() +
  theme(aspect.ratio = 1) + ylab("% male depth at male SNPs (RNA/WGS)") + xlab('') +
  scale_fill_brewer(type='qual', name='GEP')


g <- ggplot(df, aes(pct.male.dna,pct.male.tx, label=seqnames)) +
  geom_point(aes(color=gep)) +
  facet_wrap(~gep, scales = 'free') +
  stat_cor(method = 'spearman') + 
  theme_classic() +
  ggrepel::geom_text_repel(force = 2) +
  geom_abline(slope=1, intercept=0, color='lightgray', linetype='dashed') +
  theme(aspect.ratio = 1) +
  scale_color_brewer(type='qual', name='GEP') +
  scale_x_continuous(labels=scales::percent) +
  scale_y_continuous(labels=scales::percent) +
  xlab('WGS mean male SNP depth') +
  ylab('RNA-seq mean male SNP depth')  +
  guides(color=F)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

agg_png(snakemake@output[['png0']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g0)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
saveRDS(g0,snakemake@output[['ggp0']])
write_tsv(df,snakemake@output[['dat']])
