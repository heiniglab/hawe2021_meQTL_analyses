#' -----------------------------------------------------------------------------
#' Regress houseman groups from both expresssion and methylation data
#' 
#' @author Katharina Schmid
#'
#' @date 2021-01-12
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Getting parameters.")
# ------------------------------------------------------------------------------

# input files
fdata <- snakemake@input$data
fhouseman <- snakemake@input$houseman

# filtered and annotated output file
fout <- snakemake@output$data

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------

#Load files
expr<-readRDS(fdata)
covars<-read.table(fhouseman,sep="\t",header=TRUE)

#Regress for each gene/CpG
for(element in rownames(expr)){

  model<-lm(expr[element,]~.,data=covars, na.action=na.exclude)

  #Save always the residuals
  expr[element,]<-resid(model)

}
writeRDS(expr,fout)