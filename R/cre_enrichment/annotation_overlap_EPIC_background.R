#' -----------------------------------------------------------------------------
#' Sample a single background set CpGs-SNP pairs for the EPIC cis results.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Apr 27 15:59:05 2020
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)
library(GenomicRanges)
library(parallel)
source("R/lib.R")
source("R/cre_enrichment/snp_cpg_annotation_overlap_methods.R")

threads <- 6

# number of background samples to generate
sample_size <- 20

# the first parameter gives the iteration
iter <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

# ------------------------------------------------------------------------------
print("Load data.")

# based on the EUR EPIC discovery pairs
epic <- readRDS("data/current/meQTLs/epic_associations.rds")
snp_annot_epic <- readRDS("results/current/allele_frequencies/kora_epic.rds")
snp_annot_by_chr <- split(snp_annot_epic, snp_annot_epic$chr)
cpg_annot_epic <- readRDS("results/current/cpg_msd/kora_epic.rds")

# filter for cis only, adjust a bit
epic_cis <- epic %>% filter(`p-value` < 1e-14 & category == "cis") %>%
  mutate(cpg.pos = cpg.start) %>%
  dplyr::rename(snp = SNP)

# define background cpgs and snps
# we do a single background sampling to get an estimate for background overlap
# of chromHMM states
print("Defining background.")
epic_bg_cpgs <- filter(cpg_annot_epic, !id %in% epic$cpg) %>%
  mutate(pos=start) %>% as.data.frame()
rownames(epic_bg_cpgs) <- epic_bg_cpgs$id
epic_bg_snps <- filter(snp_annot_epic, !rsid %in% epic$SNP) %>% as.data.frame()
rownames(epic_bg_snps) <- epic_bg_snps$rsid

# needed in the matching method as a data frame with rownames
# change in method would have necessitated larger changes in the lib.R
cpg_annot_epic <- as.data.frame(cpg_annot_epic)
rownames(cpg_annot_epic) <- cpg_annot_epic$id

# define states
enhancer_states <- c("Genic enhancers", "Enhancers", "Bivalent Enhancer")
promoter_states <- c("Active TSS", "Flanking Active TSS")
prom_enh_states <- c(enhancer_states, promoter_states)

# ------------------------------------------------------------------------------
print("Sampling background.")

set.seed(42)

# the results from the 'observed' analysis
enrichment_observed <- 
  read_tsv("results/current/cre_enrichment/450k_vs_epic/annotation_overlap_epic.tsv")
epic_cis_filtered <- filter(epic_cis, cpg %in% enrichment_observed$cpg)
rm(epic, epic_cis)

# the observed CpGs for which we want to get a background
observed_full <- sample(unique(epic_cis_filtered$cpg))

# background CpGs should also be in prom/enh state (same as in the 'obesrved'
# analysis)
print("Getting states for background CpGs and subsetting them..")
fstates <- "results/current/cre_enrichment/450k_vs_epic/epic_bg_cpg_states.rds"
if(file.exists(fstates)) {
  epic_bg_cpg_states <- readRDS(fstates)
} else {
  epic_bg_cpgs_ranges <- with(cpg_annot_epic[epic_bg_cpgs$id, ],
                              GRanges(chr, IRanges(start, width = 2)))
  names(epic_bg_cpgs_ranges) <- epic_bg_cpgs$id
  
  # get states
  epic_bg_cpg_states <-
    get_weighted_chromhmm_annotation(epic_bg_cpgs_ranges)
  epic_bg_cpg_states <-
    apply(epic_bg_cpg_states, 1, function(x)
      names(x[which.max(x)]))
  
  # save for later use
  saveRDS(epic_bg_cpg_states, fstates)
}

# annotate and filter
epic_bg_cpgs <- mutate(epic_bg_cpgs, state = epic_bg_cpg_states[id]) %>%
  filter(state %in% prom_enh_states)
rownames(epic_bg_cpgs) <- epic_bg_cpgs$id

# subset according to iteration and sample size
observed <- observed_full[((iter-1)*sample_size+1):(iter*sample_size)]

# clean up
gc(full=T)

result <- lapply(observed, function(cpg) {
  print(paste0(cpg))
  snps <- filter(epic_cis_filtered, cpg == !!cpg) %>% dplyr::pull(snp)
  get_matched_region_epic(
    cpg,
    snps,
    epic_cis_filtered,
    snp_annot_by_chr,
    cpg_annot_epic,
    epic_bg_snps,
    epic_bg_cpgs, 
    threads = threads
  )
})

# ------------------------------------------------------------------------------
print("Sampling done. Getting state overlaps.")

# check state annotation for all sampled entities
sampled_cpgs <- unique(unlist(lapply(result, "[[", "cpg")))
sampled_cpg_ranges <- with(cpg_annot_epic[sampled_cpgs,], 
                           GRanges(chr, IRanges(start, width=2)))
names(sampled_cpg_ranges) <- sampled_cpgs

sampled_snps <- unique(unlist(lapply(result, "[[", "snps")))
snp_annot_epic_adj <- as.data.frame(snp_annot_epic)
rownames(snp_annot_epic_adj) <- snp_annot_epic_adj$rsid
sampled_snp_ranges <- with(snp_annot_epic_adj[sampled_snps,], 
                           GRanges(chr, IRanges(pos, width=1)))
names(sampled_snp_ranges) <- sampled_snps

sampled_cpg_states <- get_weighted_chromhmm_annotation(sampled_cpg_ranges)
sampled_cpg_states <- apply(sampled_cpg_states,1,function(x) names(x[which.max(x)]) )
sampled_snp_states <- get_weighted_chromhmm_annotation(sampled_snp_ranges)
sampled_snp_states <- apply(sampled_snp_states,1,function(x) names(x[which.max(x)]) )

# ------------------------------------------------------------------------------
print("Creating final results table.")

final <- lapply(result, function(r) {
  
  cpg <- r$cpg
  snps <- r$snps
  
  cpg_state <- sampled_cpg_states[cpg]
  snp_states <- sampled_snp_states[snps]
  
  # get states for CpG
  cpg_enh <- cpg_state %in% enhancer_states
  cpg_prom <- cpg_state %in% promoter_states
  
  # get states for SNPs - any must match
  enh_enh <- any(snp_states %in% enhancer_states & cpg_enh)
  ambig <- any(snp_states %in% enhancer_states & cpg_prom |
                 snp_states %in% promoter_states & cpg_enh)
  prom_prom <- any(snp_states %in% promoter_states & cpg_prom)
  
  tibble(
    bg_cpg = cpg,
    cpg_state = cpg_state,
    enh_enh = enh_enh,
    ambiguous = ambig,
    prom_prom = prom_prom,
    bg_snps = paste0(snps, collapse = ",")
  )
  
}) %>% bind_rows()

# ------------------------------------------------------------------------------
print("Saving results.")

write_tsv(final, 
          paste0("results/current/cre_enrichment/450k_vs_epic/epic_background_iter", 
          iter, ".tsv"))

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
