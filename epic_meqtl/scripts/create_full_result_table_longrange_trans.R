#' -----------------------------------------------------------------------------
#' Create the full result table from the comvined cis and trans meQTL results.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Dec  30 14:32:38 2019
#' -----------------------------------------------------------------------------


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
fannot <- snakemake@input$annot
flongrange <- snakemake@input$longrange
ftrans <- snakemake@input$trans
fsnp_pos <- snakemake@input$snp_positions
fcpg_pos <- snakemake@input$cpg_positions

# for filtering probes
fpoly_probes <- snakemake@input$polymorphic_probes
fcrosshyb_cpg <- snakemake@input$cross_hyb_cpg_probes
fcrosshyb_noncpg <- snakemake@input$cross_hyb_noncpg_probes

# outputs
fmeqtls_longrange <- snakemake@output$meqtls_longrange
fmeqtls_trans <- snakemake@output$meqtls_trans

# params
maf <- snakemake@params$post_maf
print(paste0("MAF is set to: ", maf, "."))

no_filter <- as.logical(snakemake@params$no_filter)

# ------------------------------------------------------------------------------
print("Load data.")
# ------------------------------------------------------------------------------
print("Positions...")
snp_positions <- read_tsv(fsnp_pos, col_names=F) %>%
  rename(id=X1, chr=X2, pos=X3)

cpg_positions <- read_tsv(fcpg_pos, col_names=F) %>%
  rename(id=X1, chr=X2, start=X3, stop=X4)

print("meQTLs...")
trans_meqtl <- read_tsv(ftrans)

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
combined <- trans_meqtl %>%
  filter(!cpg %in% probes_to_ignore | no_filter) %>%
  inner_join(y=snp_positions, by=c("SNP" = "id")) %>%
  inner_join(y=cpg_positions, by=c("cpg"= "id"), suffix=c(".snp", ".cpg")) %>%
  inner_join(y=annot, by = c("SNP" = "rsid")) %>%
  rename(pos.snp = pos.x, start.cpg = start, stop.cpg = stop) %>%
  select(-pos.y) %>%
  mutate(category = case_when(chr.snp == chr.cpg & abs(start.cpg - pos.snp)>1e6 ~ "longrange",
                              chr.snp != chr.cpg ~ "trans")) %>% 
  mutate(MAF = case_when(MAF > 0.5 ~ 1-MAF,
                         MAF < 0.5 ~ MAF)) %>%
  filter(MAF > maf | no_filter)

print(paste0("Total number of associations: ", nrow(combined)))
print(table(combined$category))

# ------------------------------------------------------------------------------
print("Saving individual result files.")
# ------------------------------------------------------------------------------

print("Longrange meQTLs.")
longrange <- combined %>% filter(category == "longrange")
saveRDS(longrange, fmeqtls_longrange)

print("Trans meQTLs.")
trans <- combined %>% filter(category == "trans")
saveRDS(trans, fmeqtls_trans)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
