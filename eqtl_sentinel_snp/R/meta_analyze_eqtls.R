#' ----------------------------------------------------------------------------
#' Meta analyze eQTL results from the three sub cohorts
#'
#' @author Katharina Schmid
#'
#' @date 2021-01-14
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(data.table)
library(tidyverse)
source("R/meta_analysis_methods.R")

# ------------------------------------------------------------------------------
print("Getting parameters.")

# total number of lines (w/o header) in files is: 5,738,181,750

# EUR lolipop
feur <- snakemake@input$eur
# SA lolipop
fsa <- snakemake@input$sa
# EUR kora 
feurk <- snakemake@input$eurk

# output file for meta associations
fout <- snakemake@output$result

threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Processing.")

# reading method using a filehandle and convertig to data.table
read_next <- function(handle, size) {
  # this seems to be the best option for now.
  # fread() throws an error (pot bug), and read_tsv() does not allow
  # to seek to the end of files (also reports error)
  tmp <- readLines(handle, n = size)
  fread(paste(tmp, collapse="\n"), 
        col.names = c("cpg", "gene", "beta", "tstat", "pvalue"),
        nThread = threads) %>%
    mutate(beta_se = beta / tstat)
}

feur_handle <- file(feur, open="r")
feurk_handle <- file(feurk, open="r")
fsa_handle <- file(fsa, open="r")

# read one line to avoid the header
h <- readLines(feur_handle,1)
h <- readLines(feurk_handle,1)
h <- readLines(fsa_handle,1)

# prepare main loop;
read <- T
size <- 1e7
i <- 1

ty <- Sys.time()

print("Starting loop...")
while(read) {
  
  #Show output after each line
  print(paste0("Progress (in lines): ", round(i*size, digits=2)))
  print("Time elapsed:")
  print(Sys.time() - ty)
  gc(full=T)

  # get data chunks
  print("Loading data...")
  eur <- read_next(feur_handle, size) %>%
    dplyr::rename(
      beta.eur = beta,
      beta_se.eur = beta_se,
      tstat.eur = tstat,
      pvalue.eur = pvalue
    )
  sa <- read_next(fsa_handle, size) %>%
    dplyr::rename(
      beta.sa = beta,
      beta_se.sa = beta_se,
      tstat.sa = tstat,
      pvalue.sa = pvalue
    )
  eurk <- read_next(feurk_handle, size) %>%
    dplyr::rename(
      beta.eurk = beta,
      beta_se.eurk = beta_se,
      tstat.eurk = tstat,
      pvalue.eurk = pvalue
    )
  
  # combine data.frames. no need to e.g. 'join', as all are sorted
  # by cpg/gene, i.e. read in in parallel and we can just combine.
  # keep only ids from one DF
  joined <- bind_cols(eur,
                      sa %>% dplyr::select(-cpg, -gene),
                      eurk %>% dplyr::select(-cpg, -gene))
  
  print("Getting meta results...")
  joined <- get_meta_withBeta_woFilter(joined)

  print("Saving.")
  if(i == 1) {
    fwrite(joined, fout, append = F, sep="\t",
           nThread = threads)
  } else {
    fwrite(joined, fout, append = T, sep="\t",
           nThread = threads)
  }
  
  i <- i + 1
  
  if(nrow(eur) < size) read <- F
}
print("Loop done.")
close(feur_handle)
close(feurk_handle)
close(fsa_handle)

print("Total time elapsed:")
print(Sys.time() - ty)

# ------------------------------------------------------------------------------
print("Done.SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
