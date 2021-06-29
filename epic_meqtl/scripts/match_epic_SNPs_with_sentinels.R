#' -----------------------------------------------------------------------------
#' Match EPIC pairs with the original 450k derived sentinel pairs.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Jan  28 14:32:38 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(Rsamtools)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data(Locations)
source("scripts/lib.R")

# disable progress bar
options(readr.show_progress = FALSE)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
fmeth_resid <- "results/current/methylation_residualized.rds"
fgeno_indiv <- "data/kora_genotypes/KORAS4F4_N3788_list_of_individuals.txt"
fdosage <- "results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/dosages/dosage_combined.gz"
fsentinel_pairs <- "data/benni.lehne/all_r02_110417_combined.tsv"
fepic_pairs <- "results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/all_associations_discovery_epic_specific_annotated.tsv"
fsnp_annotation <- "results/current/snp_annotation.rds"

# outputs
fresult <- "results/current/cpgs_matched_to_trans_sentinels.tsv"

# params
threads <- 10
r2_cutoff <- 0.2
dist_cutoff <- 1e6

# ------------------------------------------------------------------------------
print("Loading data..")
# ------------------------------------------------------------------------------
meth_resid <- readRDS(fmeth_resid)
geno_indiv <- read_tsv(fgeno_indiv, col_names=F) %>% pull(X1)
epic_pairs <- read_tsv(fepic_pairs, col_names = T)

# get snp ranges
snp_annotation <- readRDS(fsnp_annotation)
snp_ranges <- with(snp_annotation, GRanges(chr, IRanges(pos, width=1)))
names(snp_ranges) <- snp_annotation$rsid

sentinel_pairs <- read_tsv(fsentinel_pairs) %>% 
  filter(cpg.sentinel %in% colnames(meth_resid)) %>%
  filter(snp.sentinel %in% names(snp_ranges))

rm(snp_annotation)
gc()

# ------------------------------------------------------------------------------
print("Evaluating EPIC pairs.")
# ------------------------------------------------------------------------------

# used for TFBS enrichment; use only trans pairs
category <- "trans"

# get sentinel pair and epic TRANS subsets
sentinel_sub <- sentinel_pairs %>% filter(category == !!category)
epic_sub <- epic_pairs %>% filter(category == !!category)

# load genotypes for sentinel SNPs once
sentinel_snps <- unique(sentinel_sub %>% pull(snp.sentinel))
genotypes <- load_genotypes(fdosage, geno_indiv, snp_ranges[sentinel_snps])

# start processing, this takes a while.
result <- annotate_sentinel_snps_with_epic(sentinel_snps, genotypes, epic_sub, geno_indiv,
                                           r2_cutoff, dist_cutoff, snp_ranges, 
                                           fdosage, threads)

# get result data frame and annotated missing matches with NA
result_frame <- enframe(result, name="sentinel", value="cpgs") %>% 
  mutate(cpgs = unlist(cpgs)) %>% 
  mutate(cpgs=ifelse(cpgs=="", NA, cpgs))


# report the percentage of sentinel SNPs with a matched EPIC SNP
print("Sentinel matches:")
print(sum(!is.na(result_frame$cpgs)) / nrow(result_frame))

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
write_tsv(result_frame, fresult)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

