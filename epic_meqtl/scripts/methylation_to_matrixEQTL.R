#' -----------------------------------------------------------------------------
#' Prepare the EPIC methylation for association analysis. Loads both methylation
#' and covariate data and converts it to format suitable for MatrixEQTL
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Dec  6 14:32:38 2019
#' -----------------------------------------------------------------------------


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

data(Locations)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
fmeth <- snakemake@input$meth
fmeth_ids <- snakemake@input$meth_ids

fout_cpg_positions <- snakemake@output$cpg_pos
fout_methylation <- snakemake@output$meth

# ------------------------------------------------------------------------------
print("Load data.")
# ------------------------------------------------------------------------------
print("Methylation betas.")
load(fmeth)

meth_ids <- read_tsv(fmeth_ids, col_types=cols(X1=col_character()), 
                     col_names=F) %>% pull(X1)

# reorder according to mapping
meth <- beta_cut[,meth_ids] %>% as_tibble() %>%
  mutate(id = rownames(beta_cut)) %>%
  select(id, everything())

locs <- Locations %>% as_tibble() %>%
  dplyr::select(chr, start=pos) %>%
  mutate(id=rownames(Locations), end=start+1) %>%
  dplyr::select(id, everything())

locs <- inner_join(meth, locs, by=c("id", "id")) %>%
 select(id, chr, start, end)

print("Dimension of methylation data is (including probe ids): ")
print(dim(meth))

# ------------------------------------------------------------------------------
print("Write data.")
# ------------------------------------------------------------------------------
write_tsv(meth, fout_methylation)
write_tsv(locs, fout_cpg_positions)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

