#' -----------------------------------------------------------------------------
#' Compare the tfbs enrichment in the EPIC context
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Mar 11 13:51:40 2020
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------------------------
# load data

result_450k <- read_tsv("enrichment_chipseq_context_100_resample_10000.txt") %>%
  filter(max.q < 0.05)
result_epic_full <- read_tsv("enrichment_chipseq_epic_context_100_resample_10000_EPICFull.txt") %>%
  filter(max.q < 0.05)
result_epic_nonovel <- read_tsv("enrichment_chipseq_epic_context_100_resample_10000_EPICnoNovel.txt") %>%
  filter(max.q < 0.05)
result_450k_with_epic <- read_tsv("enrichment_chipseq_epic_context_100_resample_10000_450kwithEPIC.txt") %>%
  filter(max.q < 0.05)
result_450k_epic_bg <- read_tsv("enrichment_chipseq_epic_context_100_resample_10000_450konly.txt") %>%
  filter(max.q < 0.05)
result_450k_450k_bg <- read_tsv("enrichment_chipseq_epic_context_100_resample_10000_450kwith450kBG.txt") %>%
  filter(max.q < 0.05)

# ------------------------------------------------------------------------------
# make comparisons

compare <- function(enrichment1, enrichment2, names, bgs, comparison) {
  
  # sentinel level
  snp_inter <- length(intersect(enrichment1$sentinel, enrichment2$sentinel))
  snp_e1_spec <- length(setdiff(enrichment1$sentinel, enrichment2$sentinel))
  snp_e2_spec <- length(setdiff(enrichment2$sentinel, enrichment1$sentinel))
  
  ## TF level
  tf_inter <- length(intersect(enrichment1$tf.symbol, enrichment2$tf.symbol))
  tf_e1_spec <- length(setdiff(enrichment1$tf.symbol, enrichment2$tf.symbol))
  tf_e2_spec <- length(setdiff(enrichment2$tf.symbol, enrichment1$tf.symbol))
  
  ## pairs levels
  pairs.e1 <- unique(paste(enrichment1$sentinel, enrichment1$tf.symbol))
  pairs.e2 <- unique(paste(enrichment2$sentinel, enrichment2$tf.symbol))
  
  pairs_inter <- length(intersect(pairs.e1, pairs.e2))
  pairs_e1_spec <- length(setdiff(pairs.e1, pairs.e2))
  pairs_e2_spec <- length(setdiff(pairs.e2, pairs.e1))
  
  # print results
  
  cat("|", comparison, "|")
  cat(paste0(names[1], " / ", names[2], " | ", bgs[1], " / ", bgs[2], " | "))
  cat(snp_inter, "/", snp_e1_spec, "/", snp_e2_spec, " | ")
  cat(tf_inter, "/", tf_e1_spec, "/", tf_e2_spec, " | ")
  cat(pairs_inter, "/", pairs_e1_spec, "/", pairs_e2_spec, " | ")
  
  return(invisible(NULL))
}

# cohort effect
compare(result_450k, result_450k_450k_bg, 
        names=c("450k (cosmo)", "450k"), bgs=c("450k", "450k"),
        comparison = "cohort effect")

# background effect
compare(result_450k_450k_bg, result_450k_epic_bg, 
        names=c("450k", "450k"), bgs=c("450k", "epic"),
        comparison = "background effect")

# Effect of proxies
compare(result_450k_epic_bg, result_epic_nonovel, 
        names=c("450k", "proxies"), bgs=c("epic", "epic"), 
        comparison = "proxy loci effect")

# Effect of novel Cpg Loci
compare(result_450k_epic_bg, result_450k_with_epic, 
        names=c("450k", "450k+novel"), bgs=c("epic", "epic"),
        comparison = "novel loci effect")

# Effect of proxies and novel CpGs
compare(result_450k_epic_bg, result_epic_full, 
        names=c("450k", "EPIC"), bgs=c("epic", "epic"), 
        comparison = "proxies + novel effect")

# Full effect
compare(result_450k, result_epic_full, 
        names=c("450k (cosmo)", "EPIC"), bgs=c("450k", "epic"), 
        comparison = "full effect")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
