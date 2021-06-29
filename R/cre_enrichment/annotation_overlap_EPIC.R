#' -----------------------------------------------------------------------------
#' Perform the chromHMM annotation overlap analysis for the EPIC fine mapping
#' analysis.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Feb 26 13:31:00 2020
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)
library(GenomicRanges)
library(parallel)
source("R/lib.R")
source("R/cre_enrichment/snp_cpg_annotation_overlap_methods.R")

# ------------------------------------------------------------------------------
print("Set params.")
# ------------------------------------------------------------------------------
threads <- 6

# ------------------------------------------------------------------------------
print("Load and prepare data.")
# ------------------------------------------------------------------------------
# based on the EUR EPIC discovery pairs
epic <- readRDS("data/current/meQTLs/epic_associations.rds")
st1 <- read_tsv("data/current/meQTLs/st1_v2.txt")
snp_annot_st1 <- readRDS("results/current/allele_frequencies/kora.rds") %>%
  as_tibble()
snp_annot_epic <- readRDS("results/current/allele_frequencies/kora_epic.rds")
cpg_annot_epic <- readRDS("results/current/cpg_msd/kora_epic.rds")

# filter for cis only
st1_cis <- st1 %>% filter(p.eur.discovery < 1e-14,
                          snp.chr == cpg.chr & abs(snp.pos - cpg.pos) <= 1e6)
epic_cis <- epic %>% filter(`p-value` < 1e-14 & category == "cis")

# define background cpgs and snps
# we do a single background sampling to get an estimate for background overlap
# of chromHMM states
epic_bg_cpgs <- filter(cpg_annot_epic, !id %in% epic$cpg) %>%
  mutate(pos=start) %>% as.data.frame()
rownames(epic_bg_cpgs) <- epic_bg_cpgs$id
epic_bg_snps <- filter(snp_annot_epic, !rsid %in% epic$SNP) %>% as.data.frame()
rownames(epic_bg_snps) <- epic_bg_snps$rsid

# get cpg and snp ranges for annotation with chromHMM regions
# no 'unique' for SNPs, since there are some with same region but different id
# these will be kept as well.
snp_ranges_450k <- with(st1_cis %>% dplyr::select(snp,snp.chr,snp.pos) %>% distinct(), 
                        GRanges(paste0("chr", snp.chr), 
                                IRanges(snp.pos,width=1), id=snp))
names(snp_ranges_450k) <- snp_ranges_450k$id

cpg_ranges_450k <- unique(with(st1_cis, GRanges(paste0("chr", cpg.chr), 
                                                IRanges(cpg.pos,width=2), id=cpg)))
snp_ranges_epic <- with(epic_cis %>% dplyr::select(SNP,snp.chr,snp.pos) %>% distinct(), 
                        GRanges(snp.chr, 
                                IRanges(snp.pos,width=1), id=SNP))
names(snp_ranges_epic) <- snp_ranges_epic$id
cpg_ranges_epic <- unique(with(epic_cis, GRanges(cpg.chr, 
                                                 IRanges(cpg.start,width=2), id=cpg)))

# ------------------------------------------------------------------------------
# combine ranges and get chromHMM annotation for SNPs and CpGs
fcpg_state_annot <- "results/current/cre_enrichment/450k_vs_epic/cpg_state_annotation.rds"
if(file.exists(fcpg_state_annot)) {
  cpg_state_annotation <- readRDS(fcpg_state_annot)
} else {
  cpg_ranges <- unique(c(cpg_ranges_450k, cpg_ranges_epic))
  
  cpg_state_annotation <- get_weighted_chromhmm_annotation(cpg_ranges)
  saveRDS(cpg_state_annotation, 
          fcpg_state_annot)
  
}

fsnp_state_annot <- "results/current/cre_enrichment/450k_vs_epic/snp_state_annotation.rds"
if(file.exists(fsnp_state_annot)) {
  snp_state_annotation <- readRDS(fsnp_state_annot)
} else {
  snp_ranges <- c(snp_ranges_450k, snp_ranges_epic)
  
  snp_state_annotation <- get_weighted_chromhmm_annotation(snp_ranges)
  saveRDS(snp_state_annotation, 
          fsnp_state_annot)
}

# define states
enhancer_states <- c("Genic enhancers", "Enhancers", "Bivalent Enhancer")
promoter_states <- c("Active TSS", "Flanking Active TSS")
prom_enh_states <- c(enhancer_states, promoter_states)

# get snp and cpg state mapping and annotate
snp_states_to_st1 <- apply(snp_state_annotation[st1_cis$snp,],1,function(x) names(x[which.max(x)]) )
cpg_states_to_st1 <- apply(cpg_state_annotation[st1_cis$cpg,],1,function(x) names(x[which.max(x)]) )
snp_states_to_epic <- apply(snp_state_annotation[epic_cis$SNP,],1,function(x) names(x[which.max(x)]) )
cpg_states_to_epic <- apply(cpg_state_annotation[epic_cis$cpg,],1,function(x) names(x[which.max(x)]) )

st1_cis_filtered <- mutate(st1_cis, 
                           cpg_state = cpg_states_to_st1,
                           snp_state = snp_states_to_st1) %>% 
  filter(cpg_state %in% prom_enh_states)
epic_cis_filtered <- mutate(epic_cis, 
                            cpg_state = cpg_states_to_epic,
                            snp_state = snp_states_to_epic) %>%
  filter(cpg_state %in% prom_enh_states)

# some cleanup
rm(epic_cis, st1_cis)
gc(full=T)

# unique set of EPIC cis CpGs
epic_cis_cpgs <- unique(epic_cis_filtered$cpg)
st1_cis_cpgs <- unique(st1_cis_filtered$cpg)

# prepare helper to get results
get_annotation <- function(cpg_set, assocs,
                           enhancer_states, promoter_states,
                           threads = 1) {
  res <- mclapply(cpg_set, function(cpg) {
    
    # get associated SNPs
    subs <- assocs %>% filter(cpg == !!cpg)
    
    cpg_state <- subs %>% pull(cpg_state)
    cpg_enh <- cpg_state %in% enhancer_states
    cpg_prom <- cpg_state %in% promoter_states
    
    enh_enh <- any(subs$snp_state %in% enhancer_states & cpg_enh)
    ambig <- any(subs$snp_state %in% enhancer_states & cpg_prom |
                   subs$snp_state %in% promoter_states & cpg_enh)
    prom_prom <- any(subs$snp_state %in% promoter_states & cpg_prom)
    
    list(cpg = cpg, 
         enh_enh = (enh_enh),
         ambiguous = ambig,
         prom_prom = prom_prom)
    
  }, mc.cores=threads) %>% bind_rows
  
  # print basic information immediately
  total <- nrow(res)
  same_state <- sum(res$enh_enh, res$prom_prom)
  
  print("Fraction of pairs in the same state (enhancer/promoter)")
  print(same_state / total)
  
  return(res)
}

epic_res <- get_annotation(epic_cis_cpgs, epic_cis_filtered,
                           enhancer_states, promoter_states, threads)

st1_res <- get_annotation(st1_cis_cpgs, st1_cis_filtered,
                           enhancer_states, promoter_states, threads)

# annotation of CpG affiliation to arrays
cpg_annot_arrays <- 
  read_tsv("results/current/epic_annot/cpg_annotation_both_arrays.tsv")

# cpgs on both arrays
cpgs_both_arrays <- cpg_annot_arrays %>% filter(is_both) %>% pull(id)
# cpgs on epic array only
cpgs_epic_specific <- cpg_annot_arrays %>% filter(!is_450k) %>% pull(id)

# the correlations to check. This file is used on the LISA server with the residualized
# methylation data to get the correlations. The resulting file with the correlations
# need then to be read in here again such that we can finalize the analysis.
write_tsv(epic_cis %>% dplyr::select(cpg) %>% unique(), 
          "results/current/cre_enrichment/450k_vs_epic/epic_cis_meqtl_cpgs.txt")
write_tsv(st1_cis %>% dplyr::select(cpg) %>% unique(), 
          "results/current/cre_enrichment/450k_vs_epic/st1_cis_meqtl_cpgs.txt")

#epic_res <- read_tsv("results/current/cre_enrichment/450k_vs_epic/annotation_overlap_epic.tsv") %>%
#  select_at(1:4)
#st1_res <- read_tsv("results/current/cre_enrichment/450k_vs_epic/annotation_overlap_450k.tsv")

# read the r2 results
epic_to_st1_r2s <- read_tsv("results/current/cre_enrichment/450k_vs_epic/epic_to_st1_r2s.tsv")
epic_res_joined <- left_join(epic_res, epic_to_st1_r2s, by = c("cpg"))
epic_res_with_st1 <- epic_res_joined %>% filter(!is.na(r2))
epic_res_without_st1 <- epic_res_joined %>% filter(is.na(r2))

# join the r2 info with the overlap results and create the final table
epic_matched_same <- sum(epic_res_with_st1$prom_prom, epic_res_with_st1$enh_enh)
epic_matched_ambi <- sum(epic_res_with_st1$ambiguous)

epic_unmatched_same <- sum(epic_res_without_st1$prom_prom, epic_res_without_st1$enh_enh)
epic_unmatched_ambi <- sum(epic_res_without_st1$ambiguous)

st1_same <- sum(st1_res$prom_prom, st1_res$enh_enh)
st1_ambi <- sum(st1_res$ambiguous)

# the tables
matched_table <- matrix(c(epic_matched_same, epic_matched_ambi,
                          st1_same, st1_ambi), ncol=2, byrow = T)
unmatched_table <- matrix(c(epic_unmatched_same, epic_unmatched_ambi,
                          st1_same, st1_ambi), ncol=2, byrow = T)
rownames(matched_table) <- rownames(unmatched_table) <- c("epic","st1")
colnames(matched_table) <- colnames(unmatched_table) <- c("same state","ambiguous")

# report fisher tests
print("Results for matched EPIC CpGs:")
fisher.test(matched_table)

print("Results for unmatched EPIC CpGs:")
fisher.test(unmatched_table)

# save results
write_tsv(epic_res_joined, 
          "results/current/cre_enrichment/450k_vs_epic/annotation_overlap_epic.tsv")
write_tsv(st1_res, 
          "results/current/cre_enrichment/450k_vs_epic/annotation_overlap_450k.tsv")

# additional analysis were we add the EPIC results (in case the CpG matches to
# a ST1 CpG) to the ST1 pairs and again perofrm the overlap analysis.
st1_cis_filtered_with_r2 <- left_join(st1_cis_filtered, epic_to_st1_r2s, by=c("cpg" = "st1_cpg")) %>%
  dplyr::select(snp, cpg, snp_state, cpg_state, epic_cpg = cpg.y, epic_r2 = r2.y)
st1_cpgs <- unique(st1_cis_filtered_with_r2$cpg)

st1_res_extended <- mclapply(st1_cpgs, function(cpg) {
  
  # get associated SNPs
  subs <- st1_cis_filtered_with_r2 %>% filter(cpg == !!cpg)
  related_epic_cpgs <- subs %>% pull(epic_cpg) %>% unique
  epic_sub <- epic_cis_filtered %>% filter(cpg %in% related_epic_cpgs)
  
  cpg_state <- unique(subs %>% pull(cpg_state))
  cpg_states_epic <- unique(epic_sub %>% pull(cpg_state))
  
  # check CpG states
  cpg_enh <- cpg_state %in% enhancer_states | any(cpg_states_epic %in% enhancer_states)
  cpg_prom <- cpg_state %in% promoter_states | any(cpg_states_epic %in% promoter_states)
  
  all_snp_states <- unique(c(subs %>% pull(snp_state), 
                             epic_sub %>% pull(snp_state)))
  
  # combine with SNP state information
  enh_enh <- any(all_snp_states %in% enhancer_states) & cpg_enh
  ambig <- (any(all_snp_states %in% enhancer_states) & cpg_prom) |
                 (any(all_snp_states %in% promoter_states) & cpg_enh)
  prom_prom <- any(all_snp_states %in% promoter_states) & cpg_prom
  
  list(cpg = cpg, 
       enh_enh = enh_enh,
       ambiguous = ambig,
       prom_prom = prom_prom)
  
}, mc.cores=threads) %>% bind_rows

# print basic information immediately
total <- nrow(st1_res_extended)
same_state <- sum(st1_res_extended$enh_enh, st1_res_extended$prom_prom)

print("Fraction of pairs in the same state (enhancer/promoter)")
print((same_state) / total)

write_tsv(st1_res_extended, 
          "results/current/cre_enrichment/450k_vs_epic/annotation_overlap_450k_with_EPIC_CpGs.tsv")

# ------------------------------------------------------------------------------
# final analysis: simply check the fraction of 450k pairs in same/different state
# and the fraction of correlated epic pairs in same/different state

# get all correlated EPIC CpGs
st1_epic_all_r2s <- read_tsv("results/current/cre_enrichment/450k_vs_epic/st1_to_epic_r2s.tsv") %>%
  drop_na()
epic_cis_filtered_corr <- filter(epic_cis_filtered, cpg %in% st1_epic_all_r2s$epic_cpg)

epic_cis_cpgs_corr <- epic_cis_filtered_corr %>% pull(cpg) %>% unique
epic_res_corr <- get_annotation(epic_cis_cpgs_corr, epic_cis_filtered_corr,
                           enhancer_states, promoter_states, threads)

# get final contingency table
epic_matched_same <- sum(epic_res_corr$prom_prom | epic_res_corr$enh_enh)
epic_matched_not_same <- (nrow(epic_res_corr) - epic_matched_same)

st1_matched_same <- sum(st1_res$prom_prom | st1_res$enh_enh)
st1_matched_not_same <- nrow(st1_res) - st1_matched_same

cont <- matrix(c(st1_matched_same, st1_matched_not_same,
                 epic_matched_same, epic_matched_not_same), ncol=2, byrow = T)

colnames(cont) <- c("has SNP in same state (prom/enh)", "no SNP in same state")
rownames(cont) <- c("ST1 CpGs", "EPIC CpGs")
print(cont)
fisher.test(cont, alternative="l")

write_tsv(epic_res_corr, "results/current/cre_enrichment/450k_vs_epic/annotation_overlap_epic.tsv")

# perform a single background sampling for the epic_cis pairs and get information 
# for EPIC background

# we need to adjust DFs a bit
epic_cis_adj <- mutate(epic_cis, cpg.pos = cpg.start) %>%
  dplyr::rename(snp = SNP)
cpg_annot_epic_adj <- as.data.frame(cpg_annot_epic)
rownames(cpg_annot_epic_adj) <- cpg_annot_epic_adj$id

test_cpg <- epic_cis_adj$cpg[1]
test_snps <- filter(epic_cis_adj, cpg == test_cpg) %>% pull(snp)

snp_annot_by_chr <- split(snp_annot_epic, snp_annot_epic$chr)
test_match <- get_matched_region_epic(
  test_cpg,
  test_snps,
  epic_cis_adj,
  snp_annot_by_chr,
  cpg_annot_epic_adj,
  epic_bg_snps,
  epic_bg_cpgs
)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
