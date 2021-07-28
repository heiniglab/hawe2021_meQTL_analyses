#' -----------------------------------------------------------------------------
#' Run matrix eQTL with no consideration for cis/trans or covariates.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Dec 10 16:52:45 2019
#' -----------------------------------------------------------------------------


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(tidyverse)
library(MatrixEQTL)

# debug
print("Num Threads:")
print(RhpcBLASctl::omp_get_num_procs())
print("BLAS Num Threads:")
print(RhpcBLASctl::blas_get_num_procs())

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# output
fout_associations <- snakemake@output$associations

# input files
fdependent <- snakemake@input$dependent
findependent <- snakemake@input$independent
findependent_subset <- snakemake@input$subset

# params
threads <- snakemake@threads
pv_threshold <- 1
use_subset <- as.logical(snakemake@params$use_subset)
if(is.na(use_subset)) {
  warning("'use_subset' not specified. Using provided subset of independent entities.")
  use_subset <- T
}

keep_non_beta <- as.logical(snakemake@params$keep_non_beta)
if(is.na(keep_non_beta)) {
  warning("'keep_non_beta' not specified. Setting default to FALSE.")
  keep_non_beta <- F
}

calculate_no_fdr <- as.logical(snakemake@params$calculate_no_fdr)
if(is.na(calculate_no_fdr)) {
  warning("'calculate_no_fdr' not specified. Setting default to FALSE.")
  calculate_no_fdr <-F
}

# set openBLAS number of threads accordingly
RhpcBLASctl::blas_set_num_threads(threads)
RhpcBLASctl::omp_set_num_threads(threads)

# ------------------------------------------------------------------------------
print(paste0("Prepare sliced data: ", date()))
# ------------------------------------------------------------------------------
print("dependent data.")
dep <- SlicedData$new()
dep$fileDelimiter <- "\t"
dep$fileOmitCharacters <- "NA"
dep$fileSkipRows <- 1
dep$fileSkipColumns <- 1
dep$fileSliceSize <- 30000
dep$LoadFile(fdependent)

# we first load the data manually, subset to the necessary entities and
# convert it to asliced dataset
print("independent data.")
indep <- read_tsv(findependent, col_names=F, skip = 1)
ids <- indep %>% pull(X1)
indep <- indep %>% select(-X1) %>% data.matrix
rownames(indep) <- ids

if(use_subset) {
  samp <- read_tsv(findependent_subset, col_names=F, skip=1)
  entity_subset <- samp %>% pull(X3)
  entity_subset <- unique(c(entity_subset, samp %>% pull(X1)))
  entity_subset <- setdiff(entity_subset, NA)
  entity_subset_avail <- entity_subset[entity_subset %in% ids]
  indep <- indep[entity_subset_avail,,drop=F]
  
  if(length(entity_subset_avail) < length(entity_subset)) {
    warning("Not all entities available in data (eg CpGs with too many NAs?")
  }
}

indep_sliced <- SlicedData$new()
indep_sliced$fileOmitCharacters <- "NA" # denote missing values;
indep_sliced$fileSliceSize = 30000 # read file in pieces of 30,000 rows
indep_sliced$CreateFromMatrix(indep)

# ------------------------------------------------------------------------------
print(paste0("Compute QTLs: ", date()))
# ------------------------------------------------------------------------------

# run analysis
result <- Matrix_eQTL_main(
  snps = indep_sliced,
  gene = dep,
  output_file_name = fout_associations,
  pvOutputThreshold = pv_threshold,
  useModel = modelLINEAR,
  verbose = TRUE,
  errorCovariance = numeric(),
  pvalue.hist = FALSE,
  noFDRsaveMemory = calculate_no_fdr)

# ------------------------------------------------------------------------------
print(paste0("Analysis done: ", date()))
# ------------------------------------------------------------------------------
print("Time in seconds:")
print(result$time.in.sec)

print("Finalizing results, setting SEs.")

# if no FDR values are calculated the results are written out on disk directly and
# and not available in the output object (the following code is not working anymore)
if(! calculate_no_fdr){
  
  # set standard errors and remove unnecessary stuff if requested
  result$all$eqtls$beta_se = result$all$eqtls$beta / result$all$eqtls$statistic
  
  if(!keep_non_beta) {
    result$all$eqtls$statistic <- NULL
    result$all$eqtls$pvalue <- NULL
    result$all$eqtls$FDR <- NULL
  }

  # ------------------------------------------------------------------------------
  print("Save results.")
  # ------------------------------------------------------------------------------
  write_tsv(result$all$eqtls, path = fout_associations)
}

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

