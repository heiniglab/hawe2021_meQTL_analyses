#' -----------------------------------------------------------------------------
#' Regress covariates and houseman groups
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
fexpr <- snakemake@input$expr
fcovars <- snakemake@input$covars

# filtered and annotated output file
fout <- snakemake@output$expr

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------

#Load files
expr<-read.table(fexpr,sep="\t",header=TRUE)
covars<-read.table(fcovars,sep="\t",header=TRUE)

#Regress for each gene
for(element in rownames(expr)){

  #Perform regression for each pair
  model<-lm(expr[element,]~.,data=covars, na.action=na.exclude)

  #Save always the residuals
  expr[element,]<-resid(model)

}
write.table(expr, file=fout,
            sep="\t",quote=FALSE)