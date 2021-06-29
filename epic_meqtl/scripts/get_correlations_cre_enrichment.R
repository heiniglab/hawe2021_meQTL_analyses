library(tidyverse)
library(parallel)

source("scripts/lib.R")

# load data
epic_cpgs <- read_tsv("data/epic_cis_meqtl_cpgs.txt")
st1_cpgs <- read_tsv("data/st1_cis_meqtl_cpgs.txt")

meth <- readRDS("results/current/methylation_residualized.rds")

# get locations
epic_ranges <- get_850k_cpg_ranges()

epic_cpgs <- mutate(epic_cpgs, 
  chr = as.character(seqnames(epic_ranges[cpg])),
  pos = start(epic_ranges[cpg]))

st1_cpgs <- filter(st1_cpgs, cpg %in% names(epic_ranges)) %>%
  filter(cpg %in% colnames(meth)) %>%
  mutate(chr = as.character(seqnames(epic_ranges[cpg])),
         pos = start(epic_ranges[cpg]))


get_r2 <- function(c1, c2) {
  meth_r2 <- cor(meth[,c1], meth[,c2],
                                 use="complete")^2
  return(meth_r2)
}

threads <- 30
dist_cutoff <- 1e6
r2_cutoff <- 0.2

# we iterate separately over chromosomes
chrs <- intersect(epic_cpgs$chr, st1_cpgs$chr)

# subset methylation data
meth <- meth[,unique(c(epic_cpgs %>% pull(cpg),st1_cpgs %>% pull(cpg)))]
rm(epic_ranges)
gc(full=T)

result <- lapply(chrs, function(chr) {
  print(chr)

  # get the subset of cpgs on the current chr
  s1 <- st1_cpgs %>% filter(chr == !!chr)
  s2 <- epic_cpgs %>% filter(chr == !!chr)

  mclapply(s1$cpg, function(sc) {

    # we check only EPIC cpgs within the distance cutoff
    sc_pos <- filter(s1, cpg == !!sc) %>% pull(pos)
    s2_sub <- filter(s2, abs(pos-!!sc_pos)<=dist_cutoff) %>% pull(cpg)
#    s2_sub <- pull(s2, cpg)

    if(length(s2_sub) < 1) return(NULL)

    lapply(s2_sub, function(ec) {
      r2 <- get_r2(sc,ec)

      # if no r2 is reported, then it was too low
      if(r2 < r2_cutoff) return(NULL)

      return(list(epic_cpg = ec, st1_cpg = sc, r2=r2))

    }) %>% bind_rows # << st1 cpg lapply -----------------------------------------------

  }, mc.cores = threads) %>% bind_rows # << epic cpg mclapply --------------------------

}) %>% bind_rows # << chr lapply -------------------------------------------------------


# report results
final <- left_join(st1_cpgs, result, by=c("cpg" = "st1_cpg")) %>% dplyr::rename(st1_cpg = cpg)
write_tsv(final, "results/current/epic_st1_cre_enrichment_correlations_1MB_full.tsv")
