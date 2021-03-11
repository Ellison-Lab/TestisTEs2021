library(rtracklayer)
library(tidyverse)
library(jsonlite)
library(arrow)
library(ragg)

gtf <- import('~/work/TestisTpn/data/combined.fixed.gtf') %>%
  as_tibble()

te.lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

geps <- open_dataset('results/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

top_mods <- geps %>%
  filter(qval < optimal_ica[['qval']]) %>%
  group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
  arrange(-n_tes)

# says tep in script but have decided to show all TEs detected
tep_tes <- geps %>%
  #filter(qval < optimal_ica[['qval']]) %>%
  #filter(module == top_mods$module[1]) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  pull(gene_id) %>%
  unique()

tep_tes.no_ltrs <- tep_tes[!str_detect(tep_tes,'[-_]LTR')]

tep.name <- geps %>%
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


bws <- Sys.glob('results/finalized/bigwigs/total-rna/w1118_testes.rep*.tes.strand-forward.rpkm.bw') %>%
  set_names(.,str_extract(.,'(?<=total-rna\\/).+(?=\\.tes)')) %>%
  map(import) %>%
  map(function(x){x[seqnames(x) %in% tep_tes.no_ltrs]}) %>%
  map_df(as_tibble, .id='replicate') %>%
  group_by(seqnames, replicate) %>%
  mutate(pos = start/max(end), scaled=scale(score)[,1]) %>%
  ungroup() %>%
  mutate_if(is.factor,as.character)

g1 <- bws %>%
  ggplot(aes(pos,score, group=seqnames)) +
  geom_line(aes(color=seqnames),alpha=1) +
  theme_classic() +
  guides(color=F) +
  scale_x_continuous(name='relative position',breaks = c(0,1), labels = c('start','end')) +
  ylab('fpkm') + facet_wrap(~replicate)

x <- bws %>%
  group_by(seqnames) %>%
  mutate(pos.bin = as.factor(round(start/max(end),1.4))) %>%
  complete(replicate, seqnames, pos.bin, fill = list(score=0)) %>%
  group_by(seqnames, pos.bin) %>%
  summarise(score=mean(score,na.rm=T),.groups = 'drop_last') %>%
  mutate(scaled=scale(score)[,1])


x.hcl <- x %>% dplyr::select(seqnames, pos.bin, score) %>%
  spread(seqnames, score) %>%
  arrange(pos.bin) %>%
  as.data.frame(row.names = 'pos.bin') %>%
  column_to_rownames('pos.bin') %>% t() %>% dist() %>% hclust()

NEW.ORDER <- x.hcl$labels[x.hcl$order]

g2 <- x %>% group_by(seqnames) %>% mutate(pos.max=which.max(scaled)[1]) %>% 
  mutate(max.scaled = log2(score)) %>% ungroup() %>%
ggplot(aes(pos.bin,fct_relevel(seqnames, NEW.ORDER),fill=log2(score + 1))) +
  geom_raster(interpolate=F) +
  #scale_fill_viridis_c(name='scaled fpkm') +
  #scale_fill_distiller(type='seq',palette = 2, name='scaled fpkm', direction = 1) +
  scale_fill_gradient2(low = 'blue', mid = 'gray', high = 'red',midpoint = 0, name='log2(mean(FPKM) +1)') +
  scale_x_discrete(name='relative position',breaks = c(0,1), labels = c('start','end')) +
  ylab('')


agg_png(snakemake@output[['png1']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()


saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])
write_tsv(bws,snakemake@output[['dat']])
