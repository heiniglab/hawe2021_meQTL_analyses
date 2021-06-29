#' -----------------------------------------------------------------------------
#' Takes a dosage file and extracts positions and corresponding MAF
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jan  07 14:32:38 2019
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
fdosage_files <- snakemake@input$dosage_files

fout_annotation <- snakemake@output$annotation

# ------------------------------------------------------------------------------
print("Processing dosage files...")
# ------------------------------------------------------------------------------

all_mafs <- bind_rows(lapply(fdosage_files, function(f) {
  print(f)
  dos <- read_tsv(f, col_names=F, 
                  col_types=cols(.default=col_double(), 
                                 X1=col_character(), X2=col_character(), 
                                 X4=col_character(), X5=col_character()))
  nsamples <- (ncol(dos) - 5)
  mafs <- rowSums(select(dos, -(1:5))) / (2*nsamples)
  select(dos, c(2,1,3,4,5)) %>% 
    rename(rsid=X2, chr=X1, pos=X3, ref=X4, alt=X5) %>%
    mutate(MAF=mafs)
}))

# ------------------------------------------------------------------------------
print("Saving SNP information.")
# ------------------------------------------------------------------------------
saveRDS(all_mafs, fout_annotation)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
