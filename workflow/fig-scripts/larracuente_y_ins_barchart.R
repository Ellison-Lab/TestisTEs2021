library(tidyverse)
library(arrow)
library(ragg)
library(GenomicRanges)
library(jsonlite)
library(rtracklayer)

source("workflow/fig-scripts/theme.R")

te.lookup <- read_tsv('resources/te_id_lookup.curated.tsv.txt')

optimal_ica <- read_json('results/finalized/optimal-gep-params/larval-w1118-testes.json') %>% unlist()

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

tep.members <- w1118.gep_membership %>%
  filter(qval < optimal_ica[['qval']]) %>%
  collect() %>%
  filter(as.character(module)==tep.name) %>%
  filter(!str_detect(X1,"FBgn")) %>%
  pull(X1)

# hetchrom.ins <- open_dataset("results/finalized/hetchrom_assembly_insertions/", format='arrow')
# 
# # get 1 range for each purported insertion
# hetchrom.ins.2 <- hetchrom.ins %>%
#   collect() %>%
#   dplyr::select(chr,start,end,`repeat`,ins_id) %>%
#   distinct() %>%
#   GRanges()
# 
# hetchrom.ins.3 <- split(hetchrom.ins.2, hetchrom.ins.2$ins_id) %>%
#   lapply(FUN = function(x) GRanges(seqnames = seqnames(x), ranges=IRanges(start=min(start(x)),end = max(end(x)),names = unique(x$ins_id)),strand = strand(x)[1]) %>% .[1]) %>%
#   map_df(as_tibble, .id='ins_id') 
# 
# hetchrom.ins.4 <-  mutate(hetchrom.ins.3, `repeat`=str_extract(ins_id,regex('.+(?=\\.)'))) %>%
#   mutate(queryHits = row_number()) %>%
#   left_join(te.lookup, by=c(`repeat`='gene_id')) %>%
#   filter(is.na(component) | component %in% c(2,3,4)) %>% # removes ltr portion, because we want to test for at least potentially functional TEs
#   filter(!is.na(merged_te)) %>%
#   mutate(chrom = map_chr(as.character(seqnames), ~str_split(.,regex("_"))[[1]][1])) %>%
#   mutate(chrom = ifelse(str_detect(chrom,'^Contig'),'unmapped contig', chrom))
# 
# hetchrom.ins.5 <- hetchrom.ins.4 %>% 
#   mutate(GEP = ifelse(merged_te %in% tep.members,'TEP','other')) %>%
#   relocate(-ins_id)
# 
# df <- hetchrom.ins.5 %>%
#   mutate(chrom=ifelse(chrom=='Y','Y','other')) %>%
#   count(chrom,GEP) %>%
#   group_by(GEP) %>%
#   mutate(percent = n/sum(n))

fa <- import("resources/dmel_repbase_lib.fasta")

names(fa) <- names(fa) %>% str_remove("#.+") 

te_sizes <- width(fa) %>% set_names(names(fa)) %>% enframe(name = "te",value = "cons.length")

ins <- open_dataset("results/finalized/hetchrom_assembly_insertions/", format="arrow") %>% collect()

ins <- ins %>% left_join(te_sizes) %>% mutate(prop.missing = (cons.length - ins.size)/cons.length)

ins <- ins %>%
  left_join(te.lookup,. , by=c("gene_id"="te")) %>%
  #group_by(ins.id) %>%
  #filter(sum(del.pct) < 0.1) %>% ungroup() %>%
  mutate(GEP = ifelse(merged_te %in% tep.members,"module 27","other"))

df <- ins %>%
  filter(prop.missing < 0.1) %>%
  mutate(chrom=ifelse(str_detect(chrom,'Y'),'Y','other')) %>%
  count(chrom,GEP) %>%
  group_by(GEP) %>%
  mutate(percent = n/sum(n))

g <- df %>% 
  mutate(GEP=as.factor(GEP)) %>%
  mutate(GEP = fct_relevel(GEP,"module 27","other")) %>%
  ggplot(aes(GEP,percent,fill=chrom)) +
  geom_col(color='white') +
  theme_gte21() +
  scale_fill_gte21("binary",reverse = T) +
  theme(aspect.ratio = 1) +
  ylab('Insertion percentage') +
  scale_y_continuous(labels = scales::percent) +
  xlab("module")


has_y_ins <- ins %>% filter(str_detect(chrom,"Y") & prop.missing < 0.1) %>% pull(merged_te) %>% unique

at_least_1_df <- te.lookup %>%
  dplyr::select(merged_te) %>%
  distinct() %>%
  mutate(has.y.ins = merged_te %in% has_y_ins) %>%
  mutate(module = ifelse(merged_te %in% tep.members,"module 27","other")) 

g2 <- at_least_1_df %>%
  #filter(merged_te %in% collect(w1118.gep_membership)$X1) %>%
  group_by(module, has.y.ins) %>%
  tally() %>%
  group_by(module) %>%
  mutate(pct.has.y = n/sum(n)) %>%
  ungroup() %>%
  filter(has.y.ins) %>%
  ggplot(aes(module,pct.has.y)) +
  geom_col(color="white") +
  theme_gte21() +
  ylab('At least 1 Y insertion') +
  scale_y_continuous(labels = scales::percent)

# g2 <- ins %>%
#   filter(prop.missing < 0.1) %>%
#   dplyr::select(chrom,merged_te,ins.id,merged_te, GEP, prop.missing) %>%
#   distinct() %>%
#   left_join(te.lookup,.) %>%
#   group_by(GEP, merged_te) %>%
#   summarise(has.y.linked = any(str_detect(chrom, "Y") & prop.missing < 0.1), .groups = "drop_last") %>%
#   mutate(has.y.linked = ifelse(is.na(has.y.linked),F,has.y.linked)) %>%
#   #filter(GEP == "module 27") %>% pull(te) %>% unique()
#   #filter(GEP == "module 27") %>% #pull(te) %>% unique() %>% {tep.members[!tep.members %in% .]}
#   summarise(pct.has.y = sum(has.y.linked)/n(),n=sum(has.y.linked),n=n()) %>% 
#   ggplot(aes(GEP,pct.has.y)) +
#   geom_col(color='white') +
#   theme_gte21() +
#   ylab('At least 1 Y insertion') +
#   scale_y_continuous(labels = scales::percent) +
#   xlab("module")
  

agg_png(snakemake@output[['png']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g)
dev.off()

agg_png(snakemake@output[['png2']], width=10, height =10, units = 'in', scaling = 1.5, bitsize = 16, res = 300, background = 'transparent')
print(g2)
dev.off()

saveRDS(g,snakemake@output[['ggp']])
saveRDS(g2,snakemake@output[['ggp2']])

#write_tsv(ins,snakemake@output[['dat']])
write_tsv(at_least_1_df,snakemake@output[['dat']])


run_chisq_y <- function(top_tes,other_tes) {
  topmod_ins_gr <- top_tes[,12:14] %>% GRanges()
  othermod_ins_gr <- other_tes[,12:14] %>% GRanges()
  
  
  top_in_y<-sum(str_detect(seqnames(topmod_ins_gr),'Y'))
  
  top_not_in_y<-sum(!str_detect(seqnames(topmod_ins_gr),'Y'))
  
  other_in_y<-sum(str_detect(seqnames(othermod_ins_gr),'Y'))
  
  other_not_in_y<-sum(!str_detect(seqnames(othermod_ins_gr),'Y')) 
  
  x <- matrix(c(top_in_y, other_in_y, top_not_in_y, other_not_in_y), nrow=2,dimnames = list(c('top','other'),c('y','not_y')))
  
  list(contingency=x, tbl = broom::tidy(chisq.test(x)))
}

ins2 <- ins %>% filter(prop.missing < 0.1)

res_y <- run_chisq_y(filter(ins2,GEP=='module 27'),filter(ins2,GEP=='other'))

# Export stats info -----------------------------------------------------------------------------------

raw.stats <- res_y$tbl

stats.export <- raw.stats %>%
  mutate(script= "larracuente_y_ins_barchart.R") %>%
  mutate(desc = "overrepresentation of TEs on Y") %>%
  mutate(func = "stats::chisq.test") %>%
  mutate(ci = NA) %>%
  mutate(comparison = "TEP-TEs vs other") %>%
  mutate(alternative=NA) %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic, ci)

write_tsv(stats.export,snakemake@output[['stats']])
