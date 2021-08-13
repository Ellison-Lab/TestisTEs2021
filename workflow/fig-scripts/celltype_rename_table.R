library(tidyverse)

tibble(clusters2 = c("0/Spermatocyte",
       "1/Spermatocyte",
       "2/Spermatogonia",
       "3/Cyst",
       "4/Cyst",
       "5/Pigment",
       "6/TerminalEpithelial",
       "7/Spermatocyte",
       "8/Spermatogonia",
       "9/Spermatocyte"),
       clusters.rename = c("5/Spermatocyte","6/Spermatocyte","1/Spermatogonia",
                           "7/Cyst","8/Cyst","9/Pigment","10/TerminalEpithelial",
                           "4/Spermatocyte","2/TransitionalSpermatocyte","3/Spermatocyte")) %>%
  write_tsv(snakemake@output[['tsv']])
