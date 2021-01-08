library(tidyverse)
library(arrow)

cn <- c("gene_id", "unstranded", "first_strand", "second_strand")

files <- paste0(snakemake@input[['star']],"/ReadsPerGene.out.tab") %>%
  set_names(.,str_extract(.,regex('(?<=star\\/).+(?=\\/Reads)')))

print(files)

dat <- map(files, ~read_tsv(., col_names = cn)) %>%
	 bind_rows(.id='sample') %>%
	 filter(!gene_id %in% c("N_unmapped","N_multimapping", "N_noFeature", "N_ambiguous")) %>%
	 gather(strandedness, count, -gene_id, -sample) %>%
	 group_by(sample, strandedness)

write_dataset(dat, snakemake@output[[1]], format = "arrow")
