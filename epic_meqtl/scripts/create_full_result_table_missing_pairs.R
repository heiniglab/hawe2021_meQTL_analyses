#' -----------------------------------------------------------------------------
#' Create the full result table from the comvined cis and trans meQTL results.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Dec  30 14:32:38 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
fannot <- "results/current/snp_annotation.rds"
fmissing <- "results/current/missing_associations.tsv.gz"
fsnp_pos <- "results/current/matrixEQTL_inputs/all_snp_positions.tsv"
fcpg_pos <- "results/current/matrixEQTL_inputs/all_cpg_positions.tsv"

# for filtering probes
fpoly_probes <- "data/mccartney-et-al-2016/1-s2.0-S221359601630071X-mmc1.txt"
fcrosshyb_cpg <- "data/mccartney-et-al-2016/1-s2.0-S221359601630071X-mmc2.txt"
fcrosshyb_noncpg <- "data/mccartney-et-al-2016/1-s2.0-S221359601630071X-mmc3.txt"

# outputs
fmeqtls_combined <- "results/current/meqtls_missing_maf0.05.rds"

# params
maf <- 0.05
print(paste0("MAF is set to: ", maf, "."))

# ------------------------------------------------------------------------------
print("Load data.")
# ------------------------------------------------------------------------------
print("Positions...")
snp_positions <- read_tsv(fsnp_pos, col_names=F) %>%
  rename(id=X1, chr=X2, pos=X3)

cpg_positions <- read_tsv(fcpg_pos, col_names=F) %>%
  rename(id=X1, chr=X2, start=X3, stop=X4)

print("meQTLs...")
missing <- read_tsv(fmissing)

print("Polymorphic and cross-hybridising probes...")
polymorphic_probes <- read_tsv(fpoly_probes) %>% 
  pull(IlmnID) %>% 
  unique
print(paste0("Got ", length(polymorphic_probes), " polymorphic probes."))

cross_hyb_cpg <- read_tsv(fcrosshyb_cpg, col_names=F) %>% pull(X1) %>% unique
print(paste0("Got ", length(cross_hyb_cpg), " cross-hybridising CpG probes."))
cross_hyb_noncpg <- read_tsv(fcrosshyb_noncpg, col_names=F) %>% pull(X1) %>% unique
print(paste0("Got ", length(cross_hyb_noncpg), " cross-hybridising non-CpG probes."))

probes_to_ignore <- unique(c(polymorphic_probes, cross_hyb_cpg, cross_hyb_noncpg))
print(paste0("Number of probes to ignore: ", length(probes_to_ignore)))

annot <- readRDS(fannot)

# ------------------------------------------------------------------------------
print("Combine and add additional annotation incl. MAF.")
print("This might take a while...")
# ------------------------------------------------------------------------------
combined <- missing %>%
  filter(!cpg %in% probes_to_ignore) %>%
  inner_join(y=snp_positions, by=c("SNP" = "id")) %>%
  inner_join(y=cpg_positions, by=c("cpg"= "id"), suffix=c(".snp", ".cpg")) %>%
  inner_join(y=annot, by = c("SNP" = "rsid")) %>%
  rename(pos.snp = pos.x, start.cpg = start, stop.cpg = stop) %>%
  select(-pos.y) %>%
  mutate(category = case_when(chr.snp == chr.cpg & abs(start.cpg - pos.snp)>1e6 ~ "longrange",
                              chr.snp == chr.cpg & abs(start.cpg - pos.snp)<=1e6 ~ "cis",
                              chr.snp != chr.cpg ~ "trans")) %>% 
  mutate(MAF = case_when(MAF > 0.5 ~ 1-MAF,
                                              MAF < 0.5 ~ MAF)) %>%
  filter(MAF > maf)

# ------------------------------------------------------------------------------
print("Saving individual result files.")
# ------------------------------------------------------------------------------

print("Combined meQTLs.")
saveRDS(combined, fmeqtls_combined)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
