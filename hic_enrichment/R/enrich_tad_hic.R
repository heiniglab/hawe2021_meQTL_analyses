#' -----------------------------------------------------------------------------
#' Enrich given pairs in TADs and/or HiC data from Javierre et al.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Dec 16 16:56:46 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Loading libraries.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(GenomicRanges)
library(parallel)
library(data.table)
library(FDb.InfiniumMethylation.hg19)
source("R/lib.R")
source("R/tad_hic_enrichment/enrich_tad_hic_methods.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# inputs
fsnp_info <- snakemake@input$snp_info
fcpg_info <- snakemake@input$cpg_info
fcpgs_measured <- snakemake@input$cpgs_measured
fpairs <- snakemake@input$pairs

fmsd <- snakemake@input$msd
fallele_freq <- snakemake@input$allele_freqs
fhic <- snakemake@input$hic
fpchic <- snakemake@input$pchic
ftad <- snakemake@input$tad

# snakemake params
threads <- snakemake@threads
hic_extension <- snakemake@params$hic_ext
hic_resolution <- as.numeric(snakemake@wildcards$res)
dist_tolerance <- as.numeric(snakemake@params$dist_tolerance)
category <- snakemake@wildcards$category

# outputs
fout <- snakemake@output[[1]]

# only trans and longrange allowed
stopifnot(category %in% c("trans", "longrange"))

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------
if(category %in% "trans") {
  hic <- load_hic_regions(fhic = fhic, 
                          hic_extension, hic_resolution)
  hic_source <- hic$r1
  hic_target <- hic$r2
} else {
  hic <- load_hic_regions(fpchic = fpchic, 
                          extension = hic_extension)
  hic_source <- hic$r1
  hic_target <- hic$r2
}

if(category %in% "longrange") {
 tad_regions <- load_tad_regions(ftad)
}

# load allele frequencies and adjust colnames
snp_annotation <- readRDS(fallele_freq)
# get troubles in (global) function get.matched.snps if we have data.table
class(snp_annotation) <- "data.frame"

# load msds
cpg_annotation <- readRDS(fmsd)

# contains all SNPs/CpGs used in the analysis (also used to define background)
snp_info <- read_tsv(fsnp_info, col_names=F) %>% pull(X1)
snp_info <- intersect(snp_info, snp_annotation$rsid)
cpg_info <- read_tsv(fcpg_info, col_names=F) %>% pull(X1)
cpgs_measured <- read_tsv(fcpgs_measured, col_names=F) %>% pull(X1)

# load the actual pairs
pairs <- read_tsv(fpairs)

# ------------------------------------------------------------------------------
print("Preprocess input")
# ------------------------------------------------------------------------------
# get ranges for all SNPs
cpg_ranges <- features(FDb.InfiniumMethylation.hg19)
# for faster access later
cpg_ranges$chr <- as.character(seqnames(cpg_ranges))

snp_ranges <- with(snp_annotation, GRanges(chr, IRanges(pos, width=1)))
names(snp_ranges) <- snp_annotation$rsid

# all snps and cpgs as ranges
# create the need input list for the enrichment (list of snps and cpgs per
# sentinel row)
print("Creating input list.")
input <- mclapply(1:nrow(pairs), function(i) {
  r <- pairs[i,]
  # we only keep SNPs for which we have annotation available
  snps <- unique(c(r %>% pull("sentinel.snp"), 
                   unlist(strsplit(r %>% pull("snps"), ","))))
  snps <- intersect(snps, names(snp_ranges))
  cpgs <- unique(c(r %>% pull("sentinel.cpg"), 
                   unlist(strsplit(r %>% pull("cpgs"), ","))))
  if(length(snps) == 0 | length(cpgs) == 0) {
    NULL
  } else {
    list(snps = snps, cpgs = cpgs)
  }
}, mc.cores = threads)
# might be some SNPs were not available
input <- input[!sapply(input, is.null)]

pair_snps <- unique(unlist(lapply(input, "[[", "snps")))
pair_cpgs <- unique(unlist(lapply(input, "[[", "cpgs")))
pair_snp_ranges <- snp_ranges[pair_snps]
pair_cpg_ranges <- cpg_ranges[pair_cpgs]

# get gene annotation for promoter matching
gene_annot <- get.gene.annotation()
gene_promoters <- promoters(gene_annot)

# annotate snps and cpgs with promoter information
cpg_ranges$in_promoter <- cpg_ranges %over% gene_promoters
snp_ranges$in_promoter <- snp_ranges %over% gene_promoters


print("Creating background sets.")

# get the cpgs which we have available in our data
# DOUBLE CHECK -> do we need this the way it is here?
# i.e. it should suffice to get the ones for which we have
# msd annotation and setdiff it with the cosmo cpgs
#available_cpgs <- read.table(fcpg_ids,
#                             stringsAsFactors = F,
#                             header = F)[,1]
bg_cpgs <- setdiff(rownames(cpg_annotation), 
                   pair_cpgs)
bg_cpg_ranges <- cpg_ranges[bg_cpgs]

# snp background set
bg_snps <- setdiff(snp_annotation$rsid, snp_info)
bg_snp_ranges <- snp_ranges[bg_snps]

# ------------------------------------------------------------------------------
print("Starting calculations.")
# ------------------------------------------------------------------------------
print(paste0("Calculating overlaps using ", threads, " cores."))

print(paste0("Processing ", category, " pairs."))

if (category %in% "longrange") {
  result <- enrich_longrange(input, pair_snp_ranges, pair_cpg_ranges,
                             tad_regions,hic_source, hic_target, snp_annotation,
                             cpg_annotation, bg_cpg_ranges, bg_snp_ranges, 
                             dist_tolerance, threads)
} else if (category %in% "trans") {
  result <- enrich_trans(input, pair_snp_ranges, pair_cpg_ranges, 
                         hic_source, hic_target, 
                         snp_annotation, cpg_annotation, 
                         bg_cpg_ranges, bg_snp_ranges, threads)
}

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
save(result, file=fout)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
