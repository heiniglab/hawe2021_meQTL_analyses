#' -----------------------------------------------------------------------------
#' Regress cis eQTLs and meQTls from gene / methylation data respectively
#' 
#' @author Katharina Schmid
#'
#' @date 2021-01-12
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------

library(GenomicRanges)
source("R/scan_snps.R")

# ------------------------------------------------------------------------------
print("Getting parameters.")
# ------------------------------------------------------------------------------

# input files
fdata <- snakemake@input$data
fposition <- snakemake@input$position
fqtls <- snakemake@input$qtls
fsnps <- snakemake@input$snps
fsnps_indiv <- snakemake@input$snps_indiv

# filtered and annotated output file
fout <- snakemake@output$data

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------

#Load files
expr<-readRDS(fdata)
eqtls<-read.table(fqtls,sep="\t",header=TRUE)

#Load position information dependnt on format
if(endsWith(fposition,".RData")){
  load(fposition)
} else {
  annotation<-read.table(fposition)
}
indiv_data<-read.table(fsnps_indiv,sep="\t",header=TRUE)
  
#Cis window size
cis.window<-500000

#Iterate over each gene
for(gene in rownames(expr)){
  
  #Find correlated SNPs
  sign_snps<-eqtls$snp[eqtls$gene==gene]
  
  #Regress the values out and save the result
  if(length(sign_snps)>0){
    
    #Get all significant snps in a window of 500.000bp (cis)
    gene.annotation<-annotation[annotation$gene==gene,]
    exp.ranges<-data.frame(chr=c(paste0("chr",gene.annotation$chr)),
                           start=c(gene.annotation$start-cis.window),
                           end=c(gene.annotation$end+cis.window))
    exp.granges<-makeGRangesFromDataFrame(exp.ranges)
    
    #Extract all SNPs in the cis region of the gene
    snps<-scan.snps(exp.granges, 
                    genotype.file=fsnps, 
                    id.file=fsnps_indiv)
    
    #Find required SNP, if it is in the defined cis range
    if(sign_snps %in% rownames(snps$snps)){
      snp.data<-snps$snps[as.character(sign_snps),]
      
      #Remove individuals, where no expression was measured
      snp.data.subset<-snp.data[,colnames(snp.data) %in% indiv_data$axio_s4f4]
      
      #Reorder snp.data.subset to be in the same order as the expression values
      snp.data.subset<-snp.data.subset[,order(ordered(as.integer(
        colnames(snp.data.subset)),levels=indiv_data$axio_s4f4))]
      
      #Get gene values
      genValues<-expr[gene,]
      
      #Regress model
      model<-lm(genValues~t(snp.data.subset))
      
      #Save model results
      expr[gene,]<-model$residuals
    }
  }
}

#Save results
writeRDS(expr,fout)