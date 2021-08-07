library(rtracklayer)
library(tidyverse)
library(jsonlite)
library(arrow)
library(ragg)

source("workflow/fig-scripts/theme.R")

gtf <- import('subworkflows/gte21-custom-genome/results/custom-genome/combined.fixed.gtf') %>%
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

genes_gr <- gtf %>%
  filter(seqnames!=transcript_id) %>%
  group_by(gene_id) %>%
  filter(length(unique(transcript_id))==1) %>%
  filter(str_detect(gene_id,"FBgn") & ("mRNA" %in% type)) %>%
  filter(type %in% c("exon","5UTR","3UTR")) %>%
  split(.$gene_id) %>%
  map(~GRanges(.)) %>%
  map(~GenomicRanges::reduce(.)) %>%
  GRangesList()

tes_gr <- gtf %>% filter(gene_id %in% tes.no_ltrs) %>% GRanges() %>% split(.,.$gene_symbol) %>% GRangesList()

gr <- c(genes_gr, tes_gr)

gr <- gr %>% as.list() %>% map(~tile(.,width=50)) %>% map(unlist) %>% GRangesList() %>% unlist()

gr$gene_id <- names(gr)

bws_fw <- Sys.glob('results/finalized/bigwigs/polya-rna/larval_testes_*.strand-forward.rpkm.bw') %>%
  set_names(.,str_extract(.,'(?<=bigwigs\\/).+(?=\\.strand-)')) %>%
  map(import, which=gr) %>%
  map_df(as_tibble, .id='replicate') %>%
  mutate_if(is.factor,as.character) %>%
  mutate(strand = "-")

bws_rev <- Sys.glob('results/finalized/bigwigs/polya-rna/larval_testes_*.strand-reverse.rpkm.bw') %>%
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
  dplyr::select(replicate, chr = seqnames.x, start=start.x, end=end.x, strand=strand.x, strand.y, gene_id, score) %>%
  filter(strand==strand.y) %>%
  dplyr::select(-strand.y) %>%
  group_by(gene_id, replicate) %>%
  mutate(end = end - min(start), start = start - min(start)) %>%
  mutate(pos = start/max(end)) %>%
  ungroup() %>%
  group_by(gene_id) %>%
  filter(sum(score) > 0) %>% ungroup() %>%
  mutate(type = ifelse(str_detect(gene_id,"FBgn"),"Gene","TE"))

plot_df.summ <- plot_df %>%
  group_by(gene_id, pos, type) %>%
  summarize(score = mean(score), .groups = "drop") %>%
  mutate(bin = floor(pos * 10)/10) %>%
  group_by(gene_id, bin, type) %>%
  summarize(score=mean(score),.groups = "drop")

NEW.ORDER <- group_by(plot_df.summ, gene_id) %>% summarise(score = mean(score)) %>% arrange(score) %>% pull(gene_id)

g1 <- plot_df.summ %>%
  mutate(gene_id = fct_relevel(gene_id, NEW.ORDER)) %>%
  ggplot(aes(bin,gene_id,fill=log2(score + 1))) +
  geom_raster(interpolate=F) +
  #scale_fill_viridis_c(name='scaled fpkm') +
  scale_fill_distiller(type='seq',palette = 3, name='fpkm', direction = 1) +
  #scale_fill_gradient2(low = 'blue', mid = 'gray', high = 'red',midpoint = 0, name='log2(mean(FPKM) +1)') +
  scale_x_continuous(name='relative position', breaks = c(0,0.9), labels = c("start","end"), expand = c(0,0)) +
  ylab('') +
  facet_wrap(~type, scales = "free", strip.position="bottom") +
  theme_gte21() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank(), strip.placement = "outside")

g2 <- plot_df %>% group_by(replicate, gene_id, type) %>%
  summarise(`std. dev.` = sd(score),.groups = "drop") %>%
  ggplot(aes(type, `std. dev.`, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ggpubr::stat_compare_means(label.y = 17) +
  coord_cartesian(ylim=c(0,20)) +
  theme_gte21() +
  scale_fill_gte21("binary",reverse = T) +
  xlab("") + facet_wrap(~replicate) +
  guides(fill=F)

agg_png(snakemake@output[['png1']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g1)
dev.off()

agg_png(snakemake@output[['png2']], width=20, height =10, units = 'in', scaling = 1, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()


saveRDS(g1,snakemake@output[['ggp1']])
saveRDS(g2,snakemake@output[['ggp2']])
write_tsv(bws,snakemake@output[['dat']])


# Export stats info -----------------------------------------------------------------------------------

raw.stats <- plot_df %>% group_by(replicate, gene_id, type) %>%
  summarise(`std. dev.` = sd(score),.groups = "drop") %>%
  split(.,.$replicate) %>%
  map_df(~broom::tidy(kruskal.test(`std. dev.`~type,data=.)), .id="comparison")

stats.export <- raw.stats %>%
  mutate(script= "w1118_full_length_tx.R") %>%
  mutate(desc = "compare stdev across bins") %>%
  mutate(func = "stats::wilcox.test/ggpubr::stat_compare_means") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "positional high signal along length of TEs and genes") %>%
  mutate(alternative=NA) %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
