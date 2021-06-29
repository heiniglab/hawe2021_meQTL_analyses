#' -----------------------------------------------------------------------------
#' Match EPIC pairs with the original 450k derived sentinel pairs. Only for
#' cis meQTLs as we reinvestigate those in light of the response.
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
fmeth_resid <- snakemake@input$meth_resid
fgeno_indiv <- snakemake@input$geno_indiv
fdosage <- snakemake@input$dosage
fsentinel_pairs <- snakemake@input$sentinel_pairs
fepic_pairs <- snakemake@input$epic_pairs
fsnp_annotation <- snakemake@input$snp_annotation

# outputs
fresult <- snakemake@output$result

# params
threads <- snakemake@threads
r2_cutoff <- snakemake@params$r2_cutoff
dist_cutoff <- snakemake@params$dist_cutoff

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

# get cpg ranges
cpg_locs <- Locations %>% as_tibble() %>%
  dplyr::select(chr, start=pos) %>%
  mutate(id=rownames(Locations), end=start+1) %>%
  dplyr::select(id, everything())
cpg_ranges <- with(cpg_locs, GRanges(chr,IRanges(start,end)))
names(cpg_ranges) <- cpg_locs$id

rm(snp_annotation)
gc()

# ------------------------------------------------------------------------------
print("Evaluating EPIC pairs.")
# ------------------------------------------------------------------------------
categories <- c("longrange", "cis", "trans")
categories <- c("cis")
results <- lapply(categories, function(category) {
  print(paste0("Checking ", category, " pairs."))
  sentinel_sub <- sentinel_pairs %>% filter(category == !!category)
  epic_sub <- epic_pairs %>% filter(category == !!category)
  annotate_with_sentinels(sentinel_sub, epic_sub, 
                          meth_resid, geno_indiv,
			  r2_cutoff, dist_cutoff,
			  snp_ranges, cpg_ranges,
                          fdosage, threads)  
})

results <- bind_rows(results)

# alternative call to match only by SNPs and retrieved matched CpGs.
# used for TFBS enrichment; use only trans pairs
category <- "trans"
sentinel_sub <- sentinel_pairs %>% filter(category == !!category)
epic_sub <- epic_pairs %>% filter(category == !!category)

sentinel_snps <- unique(sentinel_sub %>% pull(snp.sentinel))
genotypes <- load_genotypes(fdosage, geno_indiv, snp_ranges[sentinel_snps])
result <- annotate_sentinel_snps_with_epic(sentinel_snps, genotypes, epic_sub, geno_indiv,
                                           r2_cutoff, dist_cutoff, snp_ranges, 
                                           fdosage, threads)

result_frame <- enframe(result, name="sentinel", value="cpgs") %>% 
  mutate(cpgs = unlist(cpgs)) %>% 
  mutate(cpgs=ifelse(cpgs=="", NA, cpgs))

write_tsv(result_frame, "results/current/cpgs_matched_to_trans_sentinels.tsv")

# report the percentage of EPIC pairs with a matched sentinel
print("Percentage of EPIC pairs with matched sentinel:")
print(sum(!is.na(results$snp.sentinel)) / nrow(results))

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(results, fresult)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

