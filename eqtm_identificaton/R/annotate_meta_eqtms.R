#' -----------------------------------------------------------------------------
#' Finalize the eQTM meta analysis. Add the eQTM category information
#' and the meta beta and sd values.
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Jul  8 08:12:09 2020
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

eqtm_input <- snakemake@input$eqtms
gene_pos_file <- snakemake@input$gene_pos
cpg_pos_file <- snakemake@input$cpg_pos

result_file <- snakemake@output$eqtms

# ------------------------------------------------------------------------------
print("Load data set and annotate it.")
# ------------------------------------------------------------------------------

# load the position information
load(cpg_pos_file) # cpg.pos dataframe
cpg_pos <- cpg.pos
gene_pos <- read_delim(gene_pos_file, delim = ",") %>%
  dplyr::select(symbol = Associated.Gene.Name,
                chr = Chromosome.Name,
                start = Gene.Start..bp.,
                end = Gene.End..bp.)

cpg_ranges <- with(cpg_pos, 
                   GRanges(chr, IRanges(start,end)))
names(cpg_ranges) <- cpg_pos$names
gene_ranges <- with(gene_pos, 
                    GRanges(paste0("chr", chr), IRanges(start, end)))
names(gene_ranges) <- gene_pos$symbol
rm(cpg_pos, gene_pos)


# load the eqtm results
gteqtm <- read_tsv(eqtm_input)

# when annotating directly the matrix eQTL results, the cpg column is called "SNP"
# correct this if necessary
if(colnames(gteqtm)[1]=="SNP"){
  print("Renaming first column from SNP to cpg!")
  colnames(gteqtm)[1]<-"cpg"
}

# get distances
dists <- distance(cpg_ranges[gteqtm$cpg], gene_ranges[gteqtm$gene]) 
gteqtm <- gteqtm %>%
  mutate(cpg_chr = as.character(seqnames(cpg_ranges[cpg])),
         cpg_pos = start(cpg_ranges[cpg]),
         gene_chr = as.character(seqnames(gene_ranges[gene])),
         gene_start = start(gene_ranges[gene]),
         gene_end = end(gene_ranges[gene]),
         category = case_when(is.na(dists) ~ "trans",
                              !is.na(dists) & abs(dists) > 1e6 ~ "longrange",
                              !is.na(dists) & abs(dists) <= 1e6 ~ "cis"))

# show the number of pairs for the individual categories
print("Nunber of pairs in categories:")
print(table(gteqtm$category))

write_tsv(gteqtm, result_file)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
