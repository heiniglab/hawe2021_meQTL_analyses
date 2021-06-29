#' -----------------------------------------------------------------------------
#' Extract common probes from 450K and EPIC arrays.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Jan  08 14:32:38 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
f450k_beta <- snakemake@input$a450k_beta
fepic_beta <- snakemake@input$epic_beta

# outputs
fout_result <- snakemake@output$result

# ------------------------------------------------------------------------------
print("Loading data..")
# ------------------------------------------------------------------------------
print("betas")
load(f450k_beta)
beta_450k <- beta
load(fepic_beta)
beta_epic <- beta

print("Getting common probes and cleaning up.")
common_probes <- intersect(rownames(beta_450k), rownames(beta_epic))

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
write_tsv(enframe(common_probes), path=fout_result, col_names=F)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
