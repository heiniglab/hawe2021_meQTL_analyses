#' -----------------------------------------------------------------------------
#' Create the full result table from all the meQTL results.
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
fmeqtls <- snakemake@input$meqtls
fsnp_pos <- snakemake@input$snp_positions
fcpg_pos <- snakemake@input$cpg_positions

# for filtering probes
fpoly_probes <- snakemake@input$polymorphic_probes
fcrosshyb_cpg <- snakemake@input$cross_hyb_cpg_probes
fcrosshyb_noncpg <- snakemake@input$cross_hyb_noncpg_probes

# outputs
fout_meqtls <- snakemake@output$meqtls

# params
maf <- snakemake@params$post_maf

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
meqtl <- read_tsv(fmeqtls)

print("Polymorphic and cross-hybridising probes...")
polymorphic_probes <- read_tsv(fpoly_probes) %>% 
  filter(EUR_AF >= 0.05) %>% 
  pull(IlmnID) %>% 
  unique
print(paste0("Got ", length(polymorphic_probes), " polymorphic probes."))

cross_hyb_cpg <- read_tsv(fcrosshyb_cpg, col_names=F) %>% pull(X1) %>% unique
print(paste0("Got ", length(cross_hyb_cpg), " cross-hybridising CpG probes."))
cross_hyb_noncpg <- read_tsv(fcrosshyb_noncpg, col_names=F) %>% pull(X1) %>% unique
print(paste0("Got ", length(cross_hyb_noncpg), " cross-hybridising non-CpG probes."))

probes_to_ignore <- unique(c(polymorphic_probes, cross_hyb_cpg, cross_hyb_noncpg))
print(paste0("Number of probes to ignore: ", length(probes_to_ignore)))

# load SNP annotation
annot <- readRDS(fannot)

# ------------------------------------------------------------------------------
print("Combine and add additional annotation incl. MAF and meQTL category.")
print("This might take a while...")
# ------------------------------------------------------------------------------
annotated <- meqtl %>%
  filter(!cpg %in% probes_to_ignore) %>%
  inner_join(y=snp_positions, by=c("SNP" = "id")) %>%
  inner_join(y=cpg_positions, by=c("cpg"= "id"), suffix=c(".snp", ".cpg")) %>%
  inner_join(y=annot, by = c("SNP" = "rsid")) %>%
  rename(snp.pos = pos.x, cpg.start = start, cpg.stop = stop) %>%
  rename(snp.chr = chr.snp, cpg.chr = chr.cpg) %>%
  select(-pos.y) %>%
#  mutate(category = "cis") %>%
  mutate(MAF = case_when(MAF > 0.5 ~ 1-MAF,
                         MAF < 0.5 ~ MAF)) %>%
  filter(MAF > maf) %>%
   # add meQTL categories
  mutate(category = case_when(snp.chr == cpg.chr & abs(cpg.start - snp.pos)>1e6 ~ "longrange",
                              snp.chr == cpg.chr & abs(cpg.start - snp.pos)<=1e6 ~ "cis",
                              snp.chr != cpg.chr ~ "trans"))


# ------------------------------------------------------------------------------
print("Saving individual result files.")
# ------------------------------------------------------------------------------
print("Total number of meQTLs:")
print(nrow(annotated))

print("Saving annotated meQTLs.")
write_tsv(annotated, fout_meqtls)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
