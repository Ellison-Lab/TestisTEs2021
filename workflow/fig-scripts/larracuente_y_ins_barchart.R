library(tidyverse)
library(arrow)
library(ragg)
library(GenomicRanges)
library(jsonlite)

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

w1118.gep_membership <- open_dataset("~/finalized/larval-w1118-testes/optimal_gep_membership", format='arrow')

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

tep.members <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  filter(as.character(module)==tep.name) %>%
  pull(X1)

hetchrom.ins <- open_dataset("~/finalized/hetchrom_assembly_insertions/", format='arrow')

# get 1 range for each purported insertion
hetchrom.ins.2 <- hetchrom.ins %>%
  collect() %>%
  dplyr::select(chr,start,end,`repeat`,ins_id) %>%
  distinct() %>%
  GRanges()

hetchrom.ins.3 <- split(hetchrom.ins.2, hetchrom.ins.2$ins_id) %>%
  lapply(FUN = function(x) GRanges(seqnames = seqnames(x), ranges=IRanges(start=min(start(x)),end = max(end(x)),names = unique(x$ins_id)),strand = strand(x)[1]) %>% .[1]) %>%
  map_df(as_tibble, .id='ins_id') 

hetchrom.ins.4 <-  mutate(hetchrom.ins.3, `repeat`=str_extract(ins_id,regex('.+(?=\\.)'))) %>%
  mutate(queryHits = row_number()) %>%
  left_join(te.lookup, by=c(`repeat`='gene_id')) %>%
  filter(is.na(component) | component %in% c(2,3,4)) %>% # removes ltr portion
  filter(!is.na(merged_te)) %>%
  mutate(chrom = map_chr(as.character(seqnames), ~str_split(.,regex("_"))[[1]][1])) %>%
  mutate(chrom = ifelse(str_detect(chrom,'^Contig'),'unmapped contig', chrom))

hetchrom.ins.5 <- hetchrom.ins.4 %>% 
  mutate(GEP = ifelse(merged_te %in% tep.members,'TEP','other')) %>%
  relocate(-ins_id)

df <- hetchrom.ins.5 %>%
  #mutate(chrom=ifelse(chrom=='Y','Y','other')) %>%
  count(chrom,GEP) %>%
  group_by(GEP) %>%
  mutate(percent = n/sum(n))

g <-ggplot(df, aes(GEP,percent,fill=chrom)) +
  geom_col(color='white') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlab('insertions') +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_viridis_d(option = 'C',direction = -1)

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
write_tsv(df,snakemake@output[['dat']])


run_chisq_y <- function(top_tes,other_tes) {
  topmod_ins_gr <- top_tes[,1:3] %>% GRanges()
  othermod_ins_gr <- other_tes[,1:3] %>% GRanges()
  
  
  top_in_y<-sum(str_detect(seqnames(topmod_ins_gr),'Y'))
  
  top_not_in_y<-sum(!str_detect(seqnames(topmod_ins_gr),'Y'))
  
  other_in_y<-sum(str_detect(seqnames(othermod_ins_gr),'Y'))
  
  other_not_in_y<-sum(!str_detect(seqnames(othermod_ins_gr),'Y')) 
  
  x <- matrix(c(top_in_y, other_in_y, top_not_in_y, other_not_in_y), nrow=2,dimnames = list(c('top','other'),c('y','not_y')))
  
  list(contingency=x, tbl = broom::tidy(chisq.test(x)))
}

res_y <- run_chisq_y(filter(hetchrom.ins.5,GEP=='TEP'),filter(hetchrom.ins.5,GEP=='other'))

