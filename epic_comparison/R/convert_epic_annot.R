#' -----------------------------------------------------------------------------
#' Annotate nearest gene distance for the EPIC based SNP and CpG annotation
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Feb 28 14:00:03 2020
#' -----------------------------------------------------------------------------
sink(file(snakemake@log[[1]], open="wt"))

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
source("R/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# input
fannot_snps <- "results/current/epic_annot/kora_epic.rds"
fannot_cpgs <- "results/current/epic_annot/cpg_annotation_both_arrays.tsv"
  
# output
fout_annot_cpgs <- "results/current/cpg_msd/kora_epic.rds"
fout_annot_snps <- "results/current/allele_frequencies/kora_epic.rds"

# ------------------------------------------------------------------------------
print("Start processing.")
# ------------------------------------------------------------------------------
snp_annot <- readRDS(fannot_snps)
# the annotation contained also 450k probes, but these do not have position 
# information annotated -> here we only want the epic probes
cpg_annot <- read_tsv(fannot_cpgs) %>% drop_na(chr)

# load gene annotaiton
ga <- get.gene.annotation()

cpg_ranges <- with(cpg_annot, GRanges(chr, IRanges(start = start, end = end)))
snp_ranges <- with(snp_annot, GRanges(chr, IRanges(pos, width=1)))

# get distance to nearest genes
distances <- distanceToNearest(snp_ranges, ga, ignore.strand = T)
distances <- as.data.frame(distances)[[3]]

snp_annot <- snp_annot %>% mutate(distGene = distances) %>%
  dplyr::rename(maf = MAF)

distances <- distanceToNearest(cpg_ranges, ga, ignore.strand = T)
distances <- as.data.frame(distances)[[3]]

cpg_annot <- cpg_annot %>% mutate(distGene = distances) %>%
  dplyr::select(id, meanBeta = mean.epic, sdBeta = sd.epic, distGene, 
                chr, start, end) %>%
  drop_na()

# ------------------------------------------------------------------------------
print("Save result.")
# ------------------------------------------------------------------------------
saveRDS(snp_annot, fout_annot_snps)
saveRDS(cpg_annot, fout_annot_cpgs)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
