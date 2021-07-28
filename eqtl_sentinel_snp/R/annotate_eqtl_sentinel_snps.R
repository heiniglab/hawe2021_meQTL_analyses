#' -----------------------------------------------------------------------------
#' Annotate eQTLs in based on distance in cis, longrange and trans
#' and filter them based on Bonferroni threshold
#' 
#' Remark: total number of tests for Bonferroni currently hardcoded!
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
library(GenomicRanges)

# ------------------------------------------------------------------------------
print("Getting parameters.")
# ------------------------------------------------------------------------------

# input files
finput <- snakemake@input$eqtls
fgenepos <- snakemake@input$gene_pos
fsnppos <- snakemake@input$snp_pos
  
# filtered and annotated output file
fout <- snakemake@output$result

threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------

#Calculate Bonferroni corrected cutoff
num_tested_pairs<- 626663229
sign.threshold <- 0.05 / num_tested_pairs

print(paste("Number of tested pairs (=number lines)",num_tested_pairs))
print("Important that this number is corrected otherwise the script is not working as supposed to!")

print(paste("Significance threshold:",sign.threshold))

#Load gene annotations and convert it to a Genomic ranges object
gene_pos <- read_delim(fgenepos, delim = ",") %>%
  dplyr::select(symbol = Associated.Gene.Name,
                chr = Chromosome.Name,
                start = Gene.Start..bp.,
                end = Gene.End..bp.)
gene_ranges <- with(gene_pos, 
                    GRanges(paste0("chr", chr), IRanges(start, end)))
names(gene_ranges) <- gene_pos$symbol
rm(gene_pos)

#Load SNP annotations
snp_pos <- read_delim(fsnppos, delim ="\t") %>%
  dplyr::select(snp.sentinel,
                snp.sentinel.chr,
                snp.sentinel.pos)
snp_pos <- unique(snp_pos)
snp_ranges <- with(snp_pos, 
                    GRanges(paste0("chr", snp.sentinel.chr), 
                            IRanges(snp.sentinel.pos, snp.sentinel.pos)))
names(snp_ranges) <- snp_pos$snp.sentinel
rm(snp_pos)

# read method using a filehandle and convertig to data.table
read_next <- function(handle, size) {
  # this seems to be the best option for now.
  # fread() throws an error (pot bug), and read_tsv() does not allow
  # to seek to the end of files (also reports error)
  tmp <- readLines(handle, n = size)
  fread(paste(tmp, collapse="\n"), 
        col.names = c("snp", "gene", "meta.beta", "meta.se", "meta.pval"),
        nThread = threads)
}

finput_handle <- file(finput, open="r")

# read one line to avoid the header
h <- readLines(finput_handle,1)

# prepare main loop;
read <- T
size <- 1e7
i <- 1
firstRows<-TRUE
  
ty <- Sys.time()

print("Starting loop...")
while(read) {

  # get data chunks
  print("Loading data...")
  eqtls <- read_next(finput_handle, size) %>%
    filter(meta.pval < sign.threshold)

  if(nrow(eqtls)>0){
    
    dist<-distance(snp_ranges[eqtls$snp], gene_ranges[eqtls$gene]) 
    eqtls <- eqtls %>%
      mutate(snp_chr = as.character(seqnames(snp_ranges[snp])),
             snp_pos = start(snp_ranges[snp]),
             gene_chr = as.character(seqnames(gene_ranges[gene])),
             gene_start = start(gene_ranges[gene]),
             gene_end = end(gene_ranges[gene]),
             category = case_when(is.na(dist) ~ "trans",
                                  !is.na(dist) & abs(dist) > 1e6 ~ "longrange",
                                  !is.na(dist) & abs(dist) <= 1e6 ~ "cis"))
    
    
    print("Saving.")
    if(firstRows) {
      fwrite(eqtls, fout, append = F, sep="\t",
             nThread = threads)
      firstRows<-FALSE
    } else {
      fwrite(eqtls, fout, append = T, sep="\t",
             nThread = threads)
    }
  
  }
  
  #Show output after each line
  print(paste0("Progress (in lines): ", round(i*size, digits=2)))
  print("Time elapsed:")
  print(Sys.time() - ty)
  gc(full=T)
  
  if(i*size >= num_tested_pairs) read <- F

  i <- i + 1
}
print("Loop done.")
close(finput_handle)

print("Total time elapsed:")
print(Sys.time() - ty)

# ------------------------------------------------------------------------------
print("Done.SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
