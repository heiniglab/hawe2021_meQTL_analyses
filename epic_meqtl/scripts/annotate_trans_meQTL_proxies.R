# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

library(tidyverse)
source("scripts/lib.R")
	
# load data
trans_epic <- read_tsv("results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/sentinel_associations_discovery_annotated.tsv")
load("results/current/methylation_residualized.rds")
meth <- readRDS("results/current/methylation_residualized.rds")
load("data/trans-cosmopairs_combined_151216.RData")

# get locations
epic_ranges <- get_850k_cpg_ranges()
epic_cpgs <- unique(trans_epic$cpg)
cosmo_cpgs <- unique(as.character(cosmo$cpg))

epic_cpgs <- tibble(cpg = epic_cpgs) %>%
  mutate(chr = as.character(seqnames(epic_ranges[cpg])),
         pos = start(epic_ranges[cpg]))

cosmo_cpgs <- tibble(cpg = cosmo_cpgs) %>%
  filter(cpg %in% names(epic_ranges)) %>%
  filter(cpg %in% colnames(meth)) %>%
  mutate(chr = as.character(seqnames(epic_ranges[cpg])),
         pos = start(epic_ranges[cpg]))

get_r2 <- function(c1, c2) {
  meth_r2 <- cor(meth[,c1], meth[,c2],
                                 use="complete")^2
  return(meth_r2)
}

# subset methylation data
meth <- meth[,unique(c(epic_cpgs %>% pull(cpg),cosmo_cpgs %>% pull(cpg)))]
rm(epic_ranges)

# we iterate separately over chromosomes
chrs <- intersect(epic_cpgs$chr, cosmo_cpgs$chr)

gc(full=T)

threads <- 50
dist_cutoff <- 1e6
r2_cutoff <- 0.2

result <- lapply(chrs, function(chr) {
  print(chr)

  # get the subset of cpgs on the current chr
  s1 <- epic_cpgs %>% filter(chr == !!chr)
  s2 <- cosmo_cpgs %>% filter(chr == !!chr)

  mclapply(s1$cpg, function(ec) {

    ec_pos <- filter(s1, cpg == !!ec) %>% pull(pos)
    s2_sub <- filter(s2, abs(pos-!!ec_pos)<=dist_cutoff) %>%
      pull(cpg)

    if(length(s2_sub) < 1) return(NULL)

    lapply(s2_sub, function(sc) {
      r2 <- get_r2(sc,ec)
      # if no r2 is reported, then it was too low
      if(r2 < r2_cutoff) return(NULL)

      return(list(epic_cpg = ec, cosmo_cpg = sc, r2=r2))

    }) %>% bind_rows # << st2 cpg lapply -----------------------------------------------
  }, mc.cores = threads) %>% bind_rows # << epic cpg mclapply --------------------------
}) %>% bind_rows()


final <- trans_epic %>% mutate(has_450k_proxy = cpg %in% result$epic_cpg)
table(final$has_450k_proxy)

write_tsv(final, "results/current/trans_meqtl_proxy_annotated.tsv")



