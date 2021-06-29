library(tidyverse)
source("scripts/lib.R")


transpairs <- as.data.frame(read_tsv("data/transpairs_r02_110117_converted_1MB.txt"), 
  stringsAsFactors=F)

# list of 115 trans hotspot sentinels
hotspot_sentinels <- read_tsv("data/sentinels.txt", col_names=F) %>% pull(X1)
# hotspot associated sentinel CpGs
sentinel_cpgs <- unique(transpairs %>%
  filter(sentinel.snp %in% hotspot_sentinels) %>% 
  pull(sentinel.cpg))

cpgs_of_interest <- read_tsv("trans_cpgs_with_tfbs_mapping.tsv", col_names=T)
mapped_sentinels <- unique(cpgs_of_interest$sentinel)

print(hotspot_sentinels[hotspot_sentinels %in% mapped_sentinels])
