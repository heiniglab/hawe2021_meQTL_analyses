#' -----------------------------------------------------------------------------
#' Load both 450k and EPIC data and create data frame with mean/sd methylation 
#' information and probe affiliation
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Jan  28 14:32:38 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(matrixStats)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# inputs
fmeth_450k <- snakemake@input$meth_450k
fmeth_epic <- snakemake@input$meth_epic

# outputs
fout_result <- snakemake@output$result

# ------------------------------------------------------------------------------
print("Load data.")
# ------------------------------------------------------------------------------
# get df with IDs and methylation mean/sd
get_df <- function(d) {
  ids <- rownames(d)
  means <- rowMeans(d, na.rm=T)
  sds <- rowSds(d, na.rm=T)

  tibble(id=ids, mean=means, sd=sds)
}

load(fmeth_450k)
df_450k <- get_df(beta)

load(fmeth_epic)
df_epic <- get_df(beta)

rm(beta)

# ------------------------------------------------------------------------------
print("Joining data and affiliating probes...")
# ------------------------------------------------------------------------------
df <- full_join(df_450k, df_epic, by = c("id" = "id"), suffix = c(".450k", ".epic")) %>%
  mutate(is_450k = !is.na(mean.450k),
         is_epic = !is.na(mean.epic),
         is_both = !is.na(mean.450k) & !is.na(mean.epic))

print("Common probes:")
print(table(df$is_both))

# annotate CpG positions for epic probes
cpg_ranges <- get_850k_cpg_ranges()
cpg_ranges_tibble <- as_tibble(as.data.frame(cpg_ranges)) %>% mutate(id=names(cpg_ranges))
df <- left_join(df, cpg_ranges_tibble %>% dplyr::select(id,chr=seqnames,start,end), by=c("id"))

# ------------------------------------------------------------------------------
print("Saving result.")
# ------------------------------------------------------------------------------
write_tsv(df, fout_result)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
