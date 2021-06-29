library(tidyverse)

fcosmo <- "data/benni.lehne/cosmopairs_combined_151216.RData"

fout_cosmo_pairs <- "results/current/cosmo_pairs.tsv"

load(fcosmo)

pairs <- cosmo %>% as_tibble() %>%
  select(snp, cpg, )

write_tsv(trans, fout_trans_pairs)
