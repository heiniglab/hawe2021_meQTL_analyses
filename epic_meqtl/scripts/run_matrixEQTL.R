#' -----------------------------------------------------------------------------
#' Run matrix eQTL on KORA genotypes and EPIC data.
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
fout_associations_cis <- snakemake@output$cis_associations
fout_associations_trans <- snakemake@output$trans_associations
fout_matrixeqtl <- snakemake@output$matrixeqtl

# input files
fmeth_data <- snakemake@input$meth
fgeno_data <- snakemake@input$geno
fcvrt <- snakemake@input$cvrt
fsnp_pos <- snakemake@input$snp_pos
fcpg_pos <- snakemake@input$cpg_pos

# params
cisdist <- as.numeric(snakemake@params$cisdist)
pvThresholdCis <- as.numeric(snakemake@params$pvThresholdCis)
pvThresholdTrans <- as.numeric(snakemake@params$pvThresholdTrans)
threads <- snakemake@threads
category <- snakemake@wildcards$category

# set openBLAS number of threads accordingly
RhpcBLASctl::blas_set_num_threads(threads)

# ------------------------------------------------------------------------------
print("Load SNP and CpG positions.")
# ------------------------------------------------------------------------------
snp_pos <- read.table(fsnp_pos, header=F, sep="\t", stringsAsFactors=F)
cpg_pos <- read.table(fcpg_pos, header=T, sep="\t", stringsAsFactors=F)

# ------------------------------------------------------------------------------
print(paste0("Prepare sliced data: ", date()))
# ------------------------------------------------------------------------------
print("Methylation data.")
# methylation data
meth <- SlicedData$new()
meth$fileDelimiter <- "\t"
meth$fileOmitCharacters <- "NA"
meth$fileSkipRows <- 1
meth$fileSkipColumns <- 1
meth$fileSliceSize <- 20000
meth$LoadFile(fmeth_data)

print("Genotype data.")
# genotype data
snps <- SlicedData$new()
snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA" # denote missing values;
snps$fileSkipRows = 0 # one row of column labels
snps$fileSkipColumns = 1 # one col of snp ids
snps$fileSliceSize = 20000 # read file in pieces of 30,000 rows
snps$LoadFile(fgeno_data)

print("Covariates.")
# covariates
cvrt <- SlicedData$new()  
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile(fcvrt)

# ------------------------------------------------------------------------------
print(paste0("Compute QTLs: ", date()))
# ------------------------------------------------------------------------------

# if a category is defined, we just touch the not needed output files to fool
# snakemake and only use the one needed
if(!is.null(category)) {
  print(paste0("Category is ", category))

  if(category %in% "cis") {
    file.create(fout_associations_trans)
    fout_associations_trans <- ""
    pvThresholdTrans <- 0
  } else {
    file.create(fout_associations_cis)
    fout_associations_cis <- ""    
    pvThresholdCis <- 0
  }  
}

# run analysis
result <- Matrix_eQTL_main(
  snps = snps,
  gene = meth,
  cvrt = cvrt,
  output_file_name.cis = fout_associations_cis,
  output_file_name = fout_associations_trans,
  pvOutputThreshold.cis = pvThresholdCis,
  pvOutputThreshold = pvThresholdTrans,
  useModel = modelLINEAR,
  snpspos = snp_pos,
  genepos = cpg_pos,
  verbose = TRUE,
  cisDist = cisdist,
  errorCovariance = numeric(),
  pvalue.hist = FALSE,
  noFDRsaveMemory = TRUE)

# ------------------------------------------------------------------------------
print(paste0("Analysis done: ", date()))
# ------------------------------------------------------------------------------
print("Time in seconds:")
print(result$time.in.sec)

# ------------------------------------------------------------------------------
print("Save results.")
# ------------------------------------------------------------------------------
saveRDS(result, fout_matrixeqtl)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

