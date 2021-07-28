#' Scans genotype files for SNPs within the provided genomic ranges
#'
#' @param ranges GRanges object containing the ranges which to scan for SNPs
#' @param genotype.file Path to and including file which contains the genotypes
#' @param id.file Path to and including file which contains the individual ids
#' 
#' @author Johann Hawe
#'
scan.snps <- function(ranges, 
                      genotype.file, 
                      id.file="individuals.txt") {
  library(data.table)

  F.SNPS = genotype.file
  F.ID = id.file
  
  # create system command using tabix...
  ranges.str <- paste(ranges, collapse=" ")
  
  cmd <- paste0("~/software/htslib-1.6/tabix ", F.SNPS, " ", ranges.str)
  data <- try(read.table(text = system(cmd,intern=TRUE) ), silent=T)
  if(inherits(data, "try-error")){
    cat("No SNPs found in specified regions.\n")
    return(list())
  } else {
    data <- data[!duplicated(data[,2]),,drop=F]
    message(paste("Processed", nrow(data), "SNPs." ))
      
    # process the genotype information to get integers/factors
    for(i in 6:ncol(data)){
      data[,i] <- factor(round(as.numeric(data[,i])))
    }
    
    ## create colnames using individual codes
    ids <- read.table(F.ID, stringsAsFactors=F, colClasses="character")
    colnames(data)<- c("chr", "name", "pos", "orig", "alt", ids[,1])
    rownames(data) <- data$name
    
    return(list(snpInfo=data[,c(1,3,4,5)], snps=data[,6:ncol(data)]))
  }
}
