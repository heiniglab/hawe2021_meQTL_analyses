#' -----------------------------------------------------------------------------
#' Run matrix eQTL with no consideration for cis/trans or covariates.
#'
#' @author Katharina Schmid
#'
#' @date Jun 4 2020
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------

library(MatrixEQTL)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

methylation_file <- snakemake@input$methylation
expression_file <- snakemake@input$expression
covariate_file <- snakemake@input$covariates
result_file <- snakemake@output$eqtms

# set openBLAS number of threads accordingly
threads <- snakemake@threads
RhpcBLASctl::blas_set_num_threads(threads)
RhpcBLASctl::omp_set_num_threads(threads)

# ------------------------------------------------------------------------------
print("Load data set.")
# ------------------------------------------------------------------------------

methylation.subset <- readRDS(methylation_file)

express.subset <- readRDS(expression_file)

expr.covars.subset <- readRDS(covariate_file)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

#Calculate all interactions (restriction to trans is not possible)
pvOutputThreshold = 1

#Load genotype data -> in our case cpgs from data matrix
snps = SlicedData$new()
snps$CreateFromMatrix(data.matrix(methylation.subset))
snps$ResliceCombined(1000)

# Load gene expression data
gene = SlicedData$new()
gene$CreateFromMatrix(data.matrix(express.subset))
gene$ResliceCombined(1000)

# Load expression covariates
cvrt = SlicedData$new()
cvrt$CreateFromMatrix(expr.covars.subset)


mytest = Matrix_eQTL_main(
  #Parameters
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = result_file,
  pvOutputThreshold    = pvOutputThreshold, 
  useModel = modelLINEAR,                          #model´s type
  errorCovariance = numeric(), #not used at the moment
  verbose = TRUE, #logical. Set to TRUE to display more detailes.
  pvalue.hist = FALSE, #"qqplot", 
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = TRUE)


# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
