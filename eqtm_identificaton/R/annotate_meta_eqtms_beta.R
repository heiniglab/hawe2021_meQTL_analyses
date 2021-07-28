#' -----------------------------------------------------------------------------
#' Add meta beta and meta se to the significant results (optional only in cis)
#'
#' @author Katharina Schmid
#'
#' @date Jul 20 2020
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------

library(tidyverse)
source("R/meta_analysis_methods.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

eqtm_input <- snakemake@input$eqtms

filter_cis <- as.logical(snakemake@params$cis)
if(length(filter_cis)==0) {
  warning("'filter_cis' not specified, setting it to FALSE.")
  filter_cis <- FALSE
}

result_file <- snakemake@output$eqtms

# ------------------------------------------------------------------------------
print("Load data set and annotate it.")
# ------------------------------------------------------------------------------

# read file
sign.eqtms<-read_tsv(eqtm_input)

# filter for cis eQTMs
if(filter_cis){
  sign.eqtms<-sign.eqtms[sign.eqtms$category=="cis",]
}

# redo meta-analysis for this data set
res<-get_meta_withBeta(sign.eqtms,filter.p=1)

# filter for required columns
res<-res[,c("cpg","gene","meta.beta","meta.sd","meta.p","category")]

write_tsv(res,result_file)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()