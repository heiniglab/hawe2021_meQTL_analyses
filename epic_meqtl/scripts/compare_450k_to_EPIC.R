#' -----------------------------------------------------------------------------
#' Compare the results obtained from the 450k and the EPIC array analyses
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Jan  08 14:32:38 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(cowplot)
library(parallel)
theme_set(theme_cowplot())

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
fcommon_probes <- snakemake@input$common_probes

f450k_disco <- snakemake@input$a450k_disco
fepic_in_st1 <- snakemake@input$epic_in_st1

# outputs
fout_missing_pairs_cis <- snakemake@output$missing_pairs_cis
fout_missing_snps_cis <- snakemake@output$missing_snps_cis
fout_missing_cpgs_cis <- snakemake@output$missing_cpgs_cis

fout_missing_pairs_longrange <- snakemake@output$missing_pairs_longrange
fout_missing_snps_longrange <- snakemake@output$missing_snps_longrange
fout_missing_cpgs_longrange <- snakemake@output$missing_cpgs_longrange

fout_missing_pairs_trans <- snakemake@output$missing_pairs_trans
fout_missing_snps_trans <- snakemake@output$missing_snps_trans
fout_missing_cpgs_trans <- snakemake@output$missing_cpgs_trans

# ------------------------------------------------------------------------------
print("Loading data..")
# ------------------------------------------------------------------------------
print("disco pairs")
disco_450k <- read_tsv(f450k_disco) %>% 
  filter(p.eur.discovery < 1e-14) %>%
  # add meQTL categories
  mutate(category = case_when(snp.chr == cpg.chr & abs(cpg.pos - snp.pos)>1e6 ~ "longrange",
                              snp.chr == cpg.chr & abs(cpg.pos - snp.pos)<=1e6 ~ "cis",
                              snp.chr != cpg.chr ~ "trans"))

print("Summary for categories in EUR disco:")
print(table(disco_450k$category))

print("epic pairs")
epic_in_st1 <- read_tsv(fepic_in_st1)

print("Getting common probes and cleaning up.")
common_probes <- read_tsv(fcommon_probes, col_names=F) %>% pull(X2)

# ------------------------------------------------------------------------------
print("Processing. Getting overlaps for common probes.")
# ------------------------------------------------------------------------------
# filter for common probes only
disco_450k_filtered <- disco_450k %>% 
  filter(cpg %in% common_probes)

epic_filtered <- epic_in_st1 %>%
  filter(cpg %in% disco_450k_filtered$cpg)

print("Joining results.")
# for each associations on the LHS, we add corresponding associations from the RHS
# if no match is possible, fields will be populated with NA for the RHS
joined <- left_join(disco_450k_filtered, epic_filtered,
                    by = c("snp" = "SNP", "cpg" = "cpg"),
                    suffix = c(".450k",".epic")) %>%
  mutate(is_epic = !is.na(beta) & `p-value` < 0.05)

print("Cleaning up a bit.")
rm(disco_450k, epic_in_st1,
   disco_450k_filtered, epic_filtered)
gc()

# get some overview numbers
total <- joined %>% tally() %>% pull(n)
total_epic <- joined %>% filter(!is.na(beta)) %>% tally() %>% pull(n)
total_cis <- joined %>% filter(category == "cis") %>% tally() %>% pull(n)
total_longrange <- joined %>% filter(category == "longrange") %>% tally() %>% pull(n)
total_trans <- joined %>% filter(category == "trans") %>% tally() %>% pull(n)

total_cis_epic <- joined %>% filter(!is.na(beta) & category == "cis") %>% tally() %>% pull(n)
total_longrange_epic <- joined %>% filter(!is.na(beta) & category == "longrange") %>% tally() %>% pull(n)
total_trans_epic <- joined %>% filter(!is.na(beta) & category == "trans") %>% tally() %>% pull(n)

print("Fractions recovered for all categories:")
print(paste0("total: ", total_epic / total, " (", total_epic, " / ", total, ")"))
print(paste0("cis: ", total_cis_epic / total_cis, " (", total_cis_epic, " / ", total_cis, ")"))
print(paste0("longrange: ", total_longrange_epic / total_longrange, " (", total_longrange_epic, " / ", total_longrange, ")"))
print(paste0("trans: ", total_trans_epic / total_trans, " (", total_trans_epic, " / ", total_trans, ")"))

# report missing pairs
missing_pairs_cis <- joined %>% filter(is.na(beta) & category == "cis")
missing_pairs_longrange <- joined %>% filter(is.na(beta) & category == "longrange")
missing_pairs_trans <- joined %>% filter(is.na(beta) & category == "trans")

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------

# cis
write_tsv(missing_pairs_cis %>% select(snp, cpg), fout_missing_pairs_cis)
write_tsv(missing_pairs_cis %>% select(snp) %>% unique(), fout_missing_snps_cis, col_names=F)
write_tsv(missing_pairs_cis %>% select(cpg) %>% unique(), fout_missing_cpgs_cis, col_names=F)

# longrange
write_tsv(missing_pairs_longrange %>% select(snp, cpg), fout_missing_pairs_longrange)
write_tsv(missing_pairs_longrange %>% select(snp) %>% unique(), fout_missing_snps_longrange, col_names=F)
write_tsv(missing_pairs_longrange %>% select(cpg) %>% unique(), fout_missing_cpgs_longrange, col_names=F)

# trans
write_tsv(missing_pairs_trans %>% select(snp, cpg), fout_missing_pairs_trans)
write_tsv(missing_pairs_trans %>% select(snp) %>% unique(), fout_missing_snps_trans, col_names=F)
write_tsv(missing_pairs_trans %>% select(cpg) %>% unique(), fout_missing_cpgs_trans, col_names=F)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
