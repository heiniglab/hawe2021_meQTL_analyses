#' -----------------------------------------------------------------------------
#' Enrich given pairs in TADs and/or HiC data from Javierre et al.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Dec 16 16:56:46 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Loading libraries.")
# ------------------------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
fpruned <- snakemake@input$pruned
fcosmo <- snakemake@input$cosmo
fout <- snakemake@output[[1]]

category <- snakemake@wildcards$category
print(paste0("category is: ", category))

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------
#load(fcosmo) # not needed just yet. we need to get the intiial mapping of
# sentinels to their original pairs.

load(fpruned)

if(category %in% "trans") {
  pruned <- trans
} else if(category %in% "longrange") {
  pruned <- longrange
} else {
  stop("category not supported.")
}
rm(trans, longrange, cis)

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------

# for now, we create a simple list with 1:1 relationships for the pruned pairs.

snps <- as.character(pruned$snp.sentinel)
cpgs <- as.character(pruned$cpg.sentinel)

to_write <- tibble(sentinel.snp = snps, snps = snps, 
                   sentinel.cpg = cpgs, cpgs = cpgs)

# ------------------------------------------------------------------------------
print("Saving pairs.")
# ------------------------------------------------------------------------------
write_tsv(to_write, fout)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
