#' -----------------------------------------------------------------------------
#' Check how many of our cis-meQTL the SNP and CpG have the same chromHMM state
#' annotation.
#' We perform random subsampling of background entities (SNPs, CpGs) to
#' create a contingency table to evaluate the overall enrichment.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Mar 15 13:35:54 2019
#' -----------------------------------------------------------------------------
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(data.table)
library(FDb.InfiniumMethylation.hg19)
library(batchtools)
source("R/lib.R")
source("R/snp_cpg_annotation_overlap_methods.R")

# ------------------------------------------------------------------------------
print("Get snakemake params")
# ------------------------------------------------------------------------------

# input
fsnp_annot <- snakemake@input$snp_annot
fcpg_annot <- snakemake@input$cpg_annot
fcosmo <- snakemake@input$cosmo

# we use these extra cosmo entity lists since the cosmo object contains
# only cis associations
fcosmo_snps <- snakemake@input$cosmo_snps
fcosmo_cpgs <- snakemake@input$cosmo_cpgs
fsnp_incl <- snakemake@input$snp_incl
fcpg_incl <- snakemake@input$cpg_incl
fcpg_msd <- snakemake@input$cpg_msd
fsnp_freq <- snakemake@input$snp_freq

# output
fout_matches <- snakemake@output$matches

# params
threads <- snakemake@threads
iteration <- as.numeric(snakemake@wildcards$iter)
samplesize <- as.numeric(snakemake@params$samplesize)
set.seed(iteration)

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------
cpg_ranges <- get450k()

# define proomter and enhancer regions
enhancer <- c("Genic enhancers", "Enhancers", "Bivalent Enhancer")
promoter <- c("Active TSS", "Flanking Active TSS")

print("Preparing cosmo object.")
load(fcosmo)
cosmo <-
  cosmo[, c("snp", "snp.chr", "snp.pos", "cpg", "cpg.chr", "cpg.pos")]
cosmo$snp <- as.character(cosmo$snp)
cosmo$cpg <- as.character(cosmo$cpg)

print("Preparing chromHMM annotations.")
env <- new.env()
load(fcpg_annot, envir = env)
cpg_annot <- with(env, weighted.annotation.all)
cpg_annot <- cpg_annot[complete.cases(cpg_annot), ]

snp_annot <- readRDS(fsnp_annot)
snp_annot <- snp_annot[complete.cases(snp_annot), ]

# get state for each CpG and SNP
cpg_states <- (sapply(1:nrow(cpg_annot), function(i) {
  return(names(which.max(cpg_annot[i, ])))
}))
names(cpg_states) <- rownames(cpg_annot)

snp_states <- unlist(sapply(1:nrow(snp_annot), function(i) {
  return(names(which.max(snp_annot[i, ])))
}))
names(snp_states) <- rownames(snp_annot)

print("Defining SNP and CpG lists.")
# get inclusion lists
load(fsnp_incl, env)
snp_incl <- with(env, snps)
load(fcpg_incl, env)
cpg_incl <- with(env, cpgs)

# load allele frequences and beta summaries
snp_freq <- readRDS(fsnp_freq)
class(snp_freq) <- "data.frame"
rownames(snp_freq) <- snp_freq$rsid
# we adapted this from a script by B. Lehne
snp_freq$maf = ifelse(snp_freq$maf > 0.5, 1 - snp_freq$maf, snp_freq$maf)
load(fcpg_msd)
cpg_msd <- msd

# only retain the cosmo pairs for which we have the freq/msd information
# available
cosmo <-
  subset(cosmo, snp %in% rownames(snp_freq) &
           cpg %in% rownames(cpg_msd))

# define background sets of CpGs and SNPs
# i.e. CpGs/SNPs not in cosmo object, but used during the analysis (meQTLs)
cosmo_snps <-
  intersect(rownames(snp_annot), fread(fcosmo_snps, header = F)$V1)
cosmo_cpgs <-
  intersect(rownames(cpg_annot), fread(fcosmo_cpgs, header = F)$V1)

cpgs_background <- setdiff(cpg_incl$rsid, cosmo_cpgs)
cpgs_background <- intersect(cpgs_background, rownames(cpg_annot))
cpgs_background <- cbind.data.frame(
  cpg = cpgs_background,
  chr = as.character(seqnames(cpg_ranges[cpgs_background])),
  pos = start(cpg_ranges[cpgs_background]),
  stringsAsFactors = F
)
rownames(cpgs_background) <- cpgs_background$cpg
cpgs_background <-
  subset(cpgs_background, cpg_states[cpg] %in% c(promoter, enhancer))

snps_background <- setdiff(snp_incl$rsid, cosmo_snps)
snps_background <- intersect(rownames(snp_annot), snps_background)
snps_background <- cbind.data.frame(
  snp = snps_background,
  chr = snp_freq[match(snps_background,
                       rownames(snp_freq)),
                 "chr"],
  pos = snp_freq[match(snps_background,
                       rownames(snp_freq)),
                 "pos"],
  stringsAsFactors = F
)
rownames(snps_background) <- snps_background$snp

print("Final preparations.")

# get all cpgs which are either in enhancer or promoter regions
cpgs <- unique(cosmo$cpg)
cpgs <-
  cbind.data.frame(cpg = cpgs,
                   state = cpg_states[match(cpgs, names(cpg_states))],
                   stringsAsFactors = F)
cpgs <- subset(cpgs, state %in% c(enhancer, promoter))

# separate by chromosome, speeds up access
snp_freq_by_chr <- split(snp_freq, snp_freq$chr)

rm(cpg_annot, snp_annot, msd, snp_freq)
gc()

# ------------------------------------------------------------------------------
print(paste0("Getting overlapping states using ", threads, " cores."))
# ------------------------------------------------------------------------------
print(paste0("Number of cpgs in enh/prom states: ", nrow(cpgs)))

nsnps_per_cpg <- tapply(cosmo$snp, cosmo$cpg, length)
qsnps <- quantile(nsnps_per_cpg, .95)
# filter out cpgs with >= 1000 snp associations (runtime issues)
snps_by_cpg <- tapply(cosmo$snp, cosmo$cpg, function(s) {
  if(length(s) <= qsnps) {
    return(s)
  } else {
    return(NULL)
  }
})
snps_by_cpg_filtered <- snps_by_cpg[!unlist(lapply(snps_by_cpg, is.null))]

cpgs_sub <- subset(cpgs, cpg %in% names(snps_by_cpg_filtered))

# sample CpGs
cpgs_to_use <- sample(cpgs_sub$cpg, samplesize)
cosmo_sub <- subset(cosmo, cpg %in% cpgs_to_use)

cosmo_sub$snp.chr <- paste0("chr", cosmo_sub$snp.chr)
cosmo_sub$cpg.chr <- paste0("chr", cosmo_sub$cpg.chr)

print(paste0("Number of cpgs being processed: ", length(cpgs_to_use)))

rm(cosmo, nsnps_per_cpg, cpgs)
gc()

cpgs_to_process <- cpgs_to_use

# take a subset for testing
#cpgs_to_process <- cpgs_to_process[sample(1:length(cpgs_to_process), 20)]

sample_background <- function(cpg_id,
                              cosmo_sub,
                              snps_by_cpg,
                              snp_freq_by_chr,
                              cpg_msd,
                              snps_background,
                              cpgs_background,
                              cpg_states,
                              snp_states,
                              enhancer,
                              promoter) {

  # get the associated snps
  cpg_snps <- snps_by_cpg[[cpg_id]]

  # begin constructing the result
  res <- c(cpg_id, paste0(cpg_snps, collapse = ","))

  # get background sample for the current CpG region
  bg_set <-
    get_matched_region(
      cpg_id,
      cpg_snps,
      cosmo_sub,
      snp_freq_by_chr,
      cpg_msd,
      snps_background,
      cpgs_background
    )

  if (is.na(bg_set) | is.null(bg_set)) {
    return(NULL)
  }

  bg_cpg <- bg_set$cpg
  bg_snps <- setdiff(bg_set$snps, NA)

  res <- c(res, bg_cpg, paste0(bg_snps, collapse = ","))

  # get annotation overlap for observed and background region
  overlap_observed <- get_overlap_region(cpg_id, cpg_snps,
                                         cpg_states, snp_states,
                                         enhancer, promoter)

  overlap_bg <-
    get_overlap_region(bg_cpg, bg_snps, cpg_states, snp_states,
                       enhancer, promoter)

  res <- c(
    res,
    overlap_observed$overlaps,
    overlap_observed$enhancer_snp,
    overlap_observed$promoter_snp,
    overlap_bg$overlaps,
    overlap_bg$enhancer_snp,
    overlap_bg$promoter_snp
  )

  return(res)
}

# define additional method arguments for running batchjobs
more.args <- list(cosmo_sub = cosmo_sub, snps_by_cpg = snps_by_cpg_filtered,
                  snp_freq_by_chr = snp_freq_by_chr, cpg_msd = cpg_msd,
                  snps_background = snps_background,
                  cpgs_background = cpgs_background,
                  cpg_states = cpg_states, snp_states = snp_states,
                  enhancer = enhancer, promoter = promoter)

res <- list(nodelist = "ibis216-010-069")

result <- run.batchtools(sample_background, cpgs_to_process, more.args = more.args,
                         name =  "state_overlap",
                         dir = paste0("results/current/chromHMM_enrichment_slurm/batchtools_", iteration),
                         n.chunks = 1000,
                         source=c("R/lib.R", "R/snp_cpg_annotation_overlap_methods.R"))

print("Get the result matrix (cpg, snp, state).")
df <- data.frame(do.call(rbind, result), stringsAsFactors = F)
colnames(df) <- c(
  "cpg", "snps",
  "bg_cpg", "bg_snps",
  "pair_enhancer",
  "pair_enh_prom",
  "pair_prom_enh",
  "pair_promoter",
  "pair_enhancer_snp",
  "pair_promoter_snp",
  "bg_enhancer",
  "bg_enh_prom",
  "bg_prom_enh",
  "bg_promoter",
  "bg_enhancer_snp",
  "bg_promoter_snp"
)

for (i in c(5:8, 11:14)) {
  df[, i] <- as.numeric(df[, i])
}

# get percentages
total <- nrow(cpgs_sub)
enh <- sum(df$pair_enhancer) / total
prom <- sum(df$pair_promoter) / total
same <- (enh + prom)
enh_prom <- sum(df$pair_enh_prom) / total
prom_enh <- sum(df$pair_prom_enh) / total
converse <- (enh_prom + prom_enh)

# print the information instantly
print("Overview over state overlap between meQTL pairs:")
print(paste0("SNP and CpG in Enhancer: ", format(enh, digits = 4)))
print(paste0("SNP and CpG in Promoter: ", format(prom, digits = 4)))
print(paste0("SNP and CpG in same annotation: ", format(same, digits = 4)))
print(paste0("Converse pairs (CpG in Enhancer): ", format(enh_prom, digits =
                                                            4)))
print(paste0("Converse pairs (CpG in Promoter): ", format(prom_enh, digits =
                                                            4)))
print(paste0("Converse pairs (total): ", format(converse, digits = 4)))

print("Saving results.")
result <- df
saveRDS(file = fout_matches,
     result)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type = "message")
