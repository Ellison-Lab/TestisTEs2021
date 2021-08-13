library(GenomicFeatures,quietly = T)
library(tidyverse, quietly = T)
library(jsonlite)
library(arrow)
library(GenomicRanges)
library(rtracklayer)

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
  pull(X1)

custom_gtf <- rtracklayer::import("subworkflows/gte21-custom-genome/results/custom-genome/combined.fixed.gtf") %>%
  as_tibble() %>%
  dplyr::select(gene_id, Blast_hit="gene_symbol") %>%
  distinct()

# ------------------------------------------
# Get IDs from larracuente supp
# ------------------------------------------

blast <- readxl::read_excel("resources/TableS8.expression.xlsx",skip = 1) %>%
  dplyr::select(gene_id_cuff="gene_id", Blast_hit, `e-value`)

preidentified <- filter(blast,str_detect(gene_id_cuff,"FBgn")) %>%
  dplyr::select(gene_id_cuff) %>%
  distinct() %>%
  left_join(custom_gtf,by=c(gene_id_cuff = "gene_id")) %>%
  mutate(gene_id = gene_id_cuff)

larrac_lookup <-filter(blast,!str_detect(gene_id_cuff,"FBgn")) %>%
  group_by(gene_id_cuff) %>%
  slice_min(`e-value`,n = 1) %>%
  dplyr::select(gene_id_cuff, Blast_hit) %>%
  distinct() %>%
  filter(n() == 1) %>%
  ungroup() %>%
  left_join(custom_gtf) %>%
  bind_rows(preidentified) %>%
  mutate(gep = ifelse(gene_id %in% tep.members,"TEP","other"))

# ------------------------------------------
# remove ambigous mappings + get introns
# ------------------------------------------

# cat Dryad.Chang.Larracuente.Genetics.2018/assembly/dmel_scaffold2_merge_V5.gtf | sed "s/ \+/\t/g" | \
#   cut -f 1,2,3,4,5,7,10,12,16 | sed 's/;//g' > chang_larracuente.tsv

gtf_df <- data.table::fread("resources/chang_larracuente.tsv",col.names = c("chr","source","type","start","end","strand","gene_id_cuff","transcript_id_cuff","weird_gene_id"),fill = T) %>%
  dplyr::select(-weird_gene_id) %>%
  as_tibble() %>%
  left_join(larrac_lookup) %>%
  mutate(gep = replace_na(gep,"unknown"))

txdb <- gtf_df %>% drop_na() %>%
  group_by(gene_id_cuff) %>%
  filter(length(unique(gene_id)) == 1) %>% # get rid of txs that map ambiguously to multiple genes
  filter(length(unique(strand)) == 1) %>% # get rid of txs that have multiple strands - this messes with the txdb function
  filter(strand %in% c("+","-")) %>% # make sure strands are known
  mutate(transcript_id = transcript_id_cuff) %>%
  GRanges() %>%  
  makeTxDbFromGRanges()

introns <- intronicParts(txdb)

y_intron_gr <- introns[str_detect(seqnames(introns),"Y")] %>% GenomicRanges::reduce(ignore.strand=T)

tep_intron_gr <- introns[as.vector(introns$gene_id) %in% tep.members] %>% GenomicRanges::reduce(ignore.strand=T)

y_tep_intron_gr <- tep_intron_gr[str_detect(seqnames(tep_intron_gr),"Y")]  %>% GenomicRanges::reduce(ignore.strand=T)

# ------------------------------------------
#  get insertions
# ------------------------------------------

# hetchrom.ins <- open_dataset("results/finalized/hetchrom_assembly_insertions/", format='arrow')
# 
# hetchrom.ins.2 <- hetchrom.ins %>%
#   collect() %>%
#   filter(str_detect(chr,"Y")) %>%
#   dplyr::select(chr,start,end,`repeat`,ins_id) %>%
#   distinct() %>%
#   GRanges()
# 
# # get 1 range for each purported insertion
# hetchrom.ins.3 <- split(hetchrom.ins.2, hetchrom.ins.2$ins_id) %>%
#   lapply(FUN = function(x) GRanges(seqnames = seqnames(x), ranges=IRanges(start=min(start(x)),end = max(end(x)),names = unique(x$ins_id)),strand = strand(x)[1]) %>% .[1]) %>%
#   map_df(as_tibble, .id='ins_id') 
# 
# hetchrom.ins.4 <-  mutate(hetchrom.ins.3, `repeat`=str_extract(ins_id,regex('.+(?=\\.)'))) %>%
#   mutate(queryHits = row_number()) %>%
#   left_join(te.lookup, by=c(`repeat`='gene_id')) %>%
#   #filter(is.na(component) | component %in% c(2,3,4)) %>% # don't filter here, need to be as loose as possible, because any solor LTR could cause TE-mapping reads
#   filter(!is.na(merged_te)) %>%
#   mutate(chrom = map_chr(as.character(seqnames), ~str_split(.,regex("_"))[[1]][1])) %>%
#   mutate(chrom = ifelse(str_detect(chrom,'^Contig'),'unmapped contig', chrom))
# 
# hetchrom.ins.5 <- hetchrom.ins.4 %>% 
#   mutate(GEP = ifelse(merged_te %in% tep.members,'TEP','other')) %>%
#   relocate(-ins_id)

fa <- import("resources/dmel_repbase_lib.fasta")

names(fa) <- names(fa) %>% str_remove("#.+") 

te_sizes <- width(fa) %>% set_names(names(fa)) %>% enframe(name = "te",value = "cons.length")

ins <- open_dataset("results/finalized/hetchrom_assembly_insertions/", format="arrow") %>% collect()

ins <- ins %>% left_join(te_sizes) %>% mutate(prop.missing = (cons.length - ins.size)/cons.length)

hetchrom.ins.5 <- ins %>%
  left_join(te.lookup, by=c(te="gene_id")) %>%
  #filter(prop.missing < 0.1) %>% # for intron fragments,even small pieces would contribute to RNA seq signal, so don't filter in this test
  mutate(GEP = ifelse(merged_te %in% tep.members,"TEP","other"))

# get a list of tes with at least 1 y linked insertion
has_y_ins <- hetchrom.ins.5 %>% filter(str_detect(chrom,"Y")) %>% pull(te) %>% unique

ins_gr <- hetchrom.ins.5 %>% 
  filter(str_detect(chrom,"Y")) %>%
  filter(te %in% has_y_ins) %>%
  dplyr::select(chrom,start,end, te, merged_te, GEP) %>% GRanges()

tep_ins <- ins_gr[ins_gr$GEP == "TEP"] %>% split(.,.$merged_te) %>% GenomicRanges::reduce(ignore.strand=T) %>% unlist()
other_ins <- ins_gr[ins_gr$GEP != "TEP"] %>% split(.,.$merged_te) %>% GenomicRanges::reduce(ignore.strand=T) %>% unlist()

tep_ins$merged_te <- names(tep_ins)
other_ins$merged_te <- names(other_ins)

# ------------------------------------------
#  make tables
# ------------------------------------------

## enrichment in any y introns

n_tep_ins_y_introns <- sum(countOverlaps(tep_ins, y_intron_gr, ignore.strand=T, type="any") >=1)
n_other_ins_y_introns <- sum(countOverlaps(other_ins, y_intron_gr, ignore.strand=T, type="any") >=1)

n_tep_ins_not_in_y_introns <- length(tep_ins) - n_tep_ins_y_introns
n_other_ins_not_in_y_introns <- length(other_ins) - n_other_ins_y_introns

ins_cont_tab <- matrix(c(n_tep_ins_y_introns, n_other_ins_y_introns, n_tep_ins_not_in_y_introns, n_other_ins_not_in_y_introns),ncol = 2)

colnames(ins_cont_tab) <- c("intron-overlapping","non-intron-overlapping")

rownames(ins_cont_tab) <- c("TEP","other")

y_intron_test <- fisher.test(ins_cont_tab)

gt_ins <- ins_cont_tab %>% 
  as_tibble(rownames = "insertions") %>%
  gt::gt() %>%
  gt::tab_header(title = "Y-linked insertions", subtitle = glue::glue("Fisher's Exact Test, OR={round(y_intron_test$estimate,digits = 2)}, P={formatC(y_intron_test$p.value,format = 'e', digits=2)}"))


## enrichment on the consensus level

n_tep_in_y_introns <- sum(countOverlaps(split(tep_ins,tep_ins$merged_te), y_intron_gr, ignore.strand=T, type="any") >= 1)
n_other_in_y_introns <- sum(countOverlaps(split(other_ins,other_ins$merged_te), y_intron_gr, ignore.strand=T, type="any") >= 1)

n_tep_not_in_y_introns <- length(unique(tep_ins$merged_te)) - n_tep_in_y_introns
n_other_not_in_y_introns <- length(unique(other_ins$merged_te)) - n_other_in_y_introns

cont_tab <- matrix(c(n_tep_in_y_introns, n_other_in_y_introns, n_tep_not_in_y_introns, n_other_not_in_y_introns),ncol = 2)

colnames(cont_tab) <- c(">0 intronic insertions","0 intronic insertions")

rownames(cont_tab) <- c("TEP","other")

y_intron_at_least_test <- fisher.test(cont_tab)

gt_consensus <- cont_tab %>% 
  as_tibble(rownames = "TE") %>%
  gt::gt() %>%
  gt::tab_header(title = "Y-linked insertions", subtitle = glue::glue("Fisher's Exact Test, OR={round(y_intron_at_least_test$estimate,digits = 2)}, P={formatC(y_intron_at_least_test$p.value,format = 'e', digits=2)}"))

## Y TEP gene intron enrichment specifically

n_tep_ins_y_tep_introns <- sum(countOverlaps(tep_ins, y_tep_intron_gr, ignore.strand=T, type="any") >=1)
n_other_ins_y_tep_introns <- sum(countOverlaps(other_ins, y_tep_intron_gr, ignore.strand=T, type="any") >=1)

n_tep_ins_not_in_y_tep_introns <- length(tep_ins) - n_tep_ins_y_tep_introns
n_other_ins_not_in_y_tep_introns <- length(other_ins) - n_other_ins_y_tep_introns

tep_ins_cont_tab <- matrix(c(n_tep_ins_y_tep_introns, n_other_ins_y_tep_introns, n_tep_ins_not_in_y_tep_introns, n_other_ins_not_in_y_tep_introns),ncol = 2)

colnames(tep_ins_cont_tab) <- c("intron-overlapping","non-intron-overlapping")

rownames(tep_ins_cont_tab) <- c("TEP","other")

y_tep_intron_test <- fisher.test(tep_ins_cont_tab)

gt_tep_ins <- tep_ins_cont_tab %>% 
  as_tibble(rownames = "insertions") %>%
  gt::gt() %>%
  gt::tab_header(title = "Y-linked insertions", subtitle = glue::glue("Fisher's Exact Test, OR={round(y_tep_intron_test$estimate,digits = 2)}, P={formatC(y_tep_intron_test$p.value,format = 'e', digits=2)}"))


gt::gtsave(gt_ins, snakemake@output[["png1"]])

gt::gtsave(gt_consensus, snakemake@output[["png2"]])

gt::gtsave(gt_tep_ins, snakemake@output[["png3"]])

# Export stats info -----------------------------------------------------------------------------------

stats.raw <- list(`enr. in chrY TEP-gene introns`=y_tep_intron_test,
     `enr. in chrY gene introns`=y_intron_test, 
     `enr. for at last 1 chrY intron insertion`=y_intron_at_least_test) %>%
  map_df(broom::tidy,.id="comparison")


stats.export <- stats.raw %>%
  mutate(script= "larracuente_y_intronic_ins.R") %>%
  mutate(desc = "enrichment in Y gene introns") %>%
  mutate(func = "stats::fisher.test") %>%
  mutate(ci = map2_chr(conf.low,conf.high,~paste(round(.x,digits = 3),round(.y,digits=3),sep=" - "))) %>%
  dplyr::select(script, comparison, desc, method, func, alternative,p.value,statistic=estimate, ci)

write_tsv(stats.export,snakemake@output[['stats']])
