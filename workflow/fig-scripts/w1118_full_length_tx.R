library(rtracklayer)
library(tidyverse)
library(jsonlite)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

gtf <- import('~/work/TestisTpn/data/combined.fixed.gtf') %>%
  as_tibble()

te.lookup <- read_tsv("resources/te_id_lookup.curated.tsv.txt")

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

geps <- open_dataset('results/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

tes <- geps %>%
  #filter(qval < optimal_ica[['qval']]) %>%
  #filter(module == top_mods$module[1]) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  filter(!str_detect(X1,'FBgn')) %>%
  pull(gene_id) %>%
  unique()

tes.no_ltrs <- tes[!str_detect(tes,'[-_]LTR')]

genes_gr <- gtf %>% filter(str_detect(gene_id,"FBgn") & type == "mRNA") %>%
  split(.$gene_id) %>%
  map(GRanges) %>%
  map(GenomicRanges::reduce) %>%
  GRangesList() %>% unlist()

tes_gr <- gtf %>% filter(gene_id %in% tes.no_ltrs) %>% GRanges()

names(tes_gr) <- tes_gr$transcript_id

gr <- c(genes_gr,tes_gr)

gr <- gr %>% tile(width=50) %>% unlist()

gr$gene_symbol <- names(gr) %>% str_extract(".+(?=\\.)")

bws_fw <- Sys.glob('~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/bigwigs/larval_testes_*.strand-forward.rpkm.bw') %>%
  set_names(.,str_extract(.,'(?<=bigwigs\\/).+(?=\\.strand-)')) %>%
  map(import, which=gr) %>%
  map_df(as_tibble, .id='replicate') %>%
  mutate_if(is.factor,as.character) %>%
  mutate(strand = "-")

bws_rev <- Sys.glob('~/amarel-scratch/TE-proj-reorg/gte21-chimeric-rnaseq/results/bigwigs/larval_testes_*.strand-reverse.rpkm.bw') %>%
  set_names(.,str_extract(.,'(?<=bigwigs\\/).+(?=\\.strand-)')) %>%
  map(import, which=gr) %>%
  map_df(as_tibble, .id='replicate') %>%
  mutate_if(is.factor,as.character) %>%
  mutate(strand="+")

bws <- bind_rows(fw = bws_fw, rev = bws_rev) %>% mutate(subjectHits = row_number())

bws_gr <- GRanges(bws)

ol_df <- findOverlaps(gr, bws_gr) %>% as_tibble()

plot_df <- as_tibble(gr) %>% 
  mutate(queryHits = row_number()) %>% 
  left_join(ol_df) %>% 
  left_join(bws, by = "subjectHits") %>% 
  dplyr::select(replicate, chr = seqnames.x, start=start.x, end=end.x, strand=strand.x, strand.y, gene_symbol, score) %>%
  filter(strand==strand.y) %>%
  dplyr::select(-strand.y) %>%
  group_by(gene_symbol, replicate) %>%
  mutate(end = end - min(start), start = start - min(start)) %>% 
  mutate(pos = start/max(end)) %>%
  ungroup() %>% 
  group_by(gene_symbol) %>%
  filter(sum(score) > 0) %>% ungroup() %>%
  mutate(type = ifelse(str_detect(gene_symbol,"FBgn"),"Gene","TE"))

plot_df.summ <- plot_df %>% 
  group_by(gene_symbol, pos, type) %>%
  summarize(score = mean(score), .groups = "drop") %>%
  mutate(bin = floor(pos * 10)/10) %>%
  group_by(gene_symbol, bin, type) %>%
  summarize(score=mean(score),.groups = "drop")

NEW.ORDER <- group_by(plot_df.summ, gene_symbol) %>% summarise(score = mean(score)) %>% arrange(score) %>% pull(gene_symbol)

g1 <- ggplot(plot_df.summ, aes(bin,fct_relevel(gene_symbol, NEW.ORDER),fill=log2(score + 1))) +
  geom_raster(interpolate=F) +
  #scale_fill_viridis_c(name='scaled fpkm') +
  scale_fill_distiller(type='seq',palette = 3, name='fpkm', direction = 1) +
  #scale_fill_gradient2(low = 'blue', mid = 'gray', high = 'red',midpoint = 0, name='log2(mean(FPKM) +1)') +
  scale_x_continuous(name='relative position', breaks = c(0,0.9), labels = c("start","end"), expand = c(0,0)) +
  ylab('') +
  facet_wrap(~type, scales = "free") +
  theme_gte21() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank())

g2 <- plot_df %>% group_by(replicate, gene_symbol, type) %>%
  summarise(`std. dev.` = sd(score),.groups = "drop") %>%
  ggplot(aes(type, `std. dev.`, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 17) +
  coord_cartesian(ylim=c(0,20)) +
  theme_gte21() +
  scale_fill_gte21("binary",reverse = T) +
  xlab("") + facet_wrap(~replicate)

agg_png(snakemake@output[['png1']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()


saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])
write_tsv(bws,snakemake@output[['dat']])
