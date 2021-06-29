#' -----------------------------------------------------------------------------
#' Enrich our identified TF list over all loci in the Lemire 2017 epigenetic
#' regulator gene list
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Apr 17 15:22:19 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(data.table)
library(ggpubr)
library(ggplot2)
library(cowplot)
source("../R/lib.R")

# ------------------------------------------------------------------------------
print("Getting snakemake info and loading data.")
# ------------------------------------------------------------------------------
fout_plot <- snakemake@output$plot
fout_data <- snakemake@output$data

flemire_regulators <- snakemake@input$lemire_regulators
flemire_tfs <- snakemake@input$lemire_tfs
fincl_snps <- snakemake@input$incl_snps
fcosmo_snps <- snakemake@input$cosmo_snps
fpruned <- snakemake@input$pruned

threads <- snakemake@threads

# load list of lemire regulators (epigenetic regulators)
lemire_regs <- unique(fread(flemire_regulators, header=F)$V1)
print("Genes in regulator set:")
print(length(lemire_regs))

lemire_tfs <- unique(fread(flemire_tfs, header=F)$V1)
print("Genes in TF set:")
print(length(lemire_tfs))

lemire_znfs <- lemire_tfs[grepl("ZNF", lemire_tfs)]
print("Genes in ZNF set:")
print(length(lemire_znfs))

# ------------------------------------------------------------------------------
print("Get gene annotation.")
# ------------------------------------------------------------------------------
ga <- get.gene.annotation()

# ------------------------------------------------------------------------------
print("Define and annotate sentinel regions.")
# ------------------------------------------------------------------------------
pairs <- fread(fpruned)
trans_pairs <- subset(pairs, chr.snp != chr.cpg)
sentinel_regions <- with(trans_pairs, GRanges(paste0("chr", chr.snp),
                                              IRanges(interval.start.snp,
                                                      interval.end.snp)))
names(sentinel_regions) <- trans_pairs$sentinel.snp
sentinel_regions <- unique(sentinel_regions)
background_sizes <- width(sentinel_regions)

obs_genes <- unlist(subsetByOverlaps(ga,sentinel_regions)$SYMBOL)
print("Number of observed genes:")
print(length(obs_genes))

# needed for background definition
load(fincl_snps) # contains 'snps' data.frame
cosmo_snps <- fread(fcosmo_snps, header=F)$V1

enrich_gene_set <- function(observed, lemire_gene_set, title="") {

  obs_lemire <- intersect(observed, lemire_gene_set)
  obs_not_lemire <- setdiff(obs_genes, lemire_gene_set)
  print("Number of observed genes in lemire set:")
  print(length(obs_lemire))
  obs_or <- length(obs_lemire) / length(obs_not_lemire)

  # ----------------------------------------------------------------------------
  print("Defining/sampling 1000 background sets and getting overlaps.")
  # ----------------------------------------------------------------------------
  # get background regions, recycle widths from sentinel regions
  overlaps <- mclapply(1:1000, function(r) {
    set.seed(r)
    background_snps <- sample(setdiff(snps$rsid, cosmo_snps),
                              size=length(sentinel_regions),
                              replace = F)
    background_snps <- snps[match(background_snps, snps$rsid),]
    background_regions <- with(background_snps, GRanges(paste0("chr", chr),
                                                      IRanges(pos, width=1)))
    background_regions <- resize(background_regions, background_sizes, fix="center")
    bg_genes <- unlist(subsetByOverlaps(ga,background_regions)$SYMBOL)
    bg_lemire <- intersect(bg_genes, lemire_gene_set)
    bg_not_lemire <- setdiff(bg_genes, lemire_gene_set)
    c(length(bg_lemire), length(bg_not_lemire))
  }, mc.cores=threads)

  bg_in_lemire <- sapply(overlaps, "[[", 1)
  bg_not_in_lemire <- sapply(overlaps, "[[", 2)

  # get the overlap ORs
  bg_ors <- bg_in_lemire / bg_not_in_lemire

  # get mean ratio for fold enrichment
  mean_bg_ors <- mean(bg_ors)

  # ----------------------------------------------------------------------------
  print("Getting ORs and plotting.")
  # ----------------------------------------------------------------------------
  ORs <- obs_or / bg_ors
  print("Fraction of background ORs leq 1 (observed over bg):")
  print((sum(ORs<=1) + 1) / (length(bg_ors)+1))
  print("Mean fold-enrichment:")
  print(obs_or / mean_bg_ors)

  # create the data frames
  toplot <- cbind.data.frame(bg_ratios=bg_ors)
  toplot2 <- cbind.data.frame(observed=obs_or)

  gp <- ggplot(toplot, aes(x=bg_ratios)) +
    geom_histogram(aes(col = "background", fill = "background"),
                   bins = length(bg_ors)/10) +
    geom_segment(
      data = toplot2,
      show.legend = F,
      aes(
        x = observed,
        xend = observed,
        y = 5,
        yend = 0,
        col = "observed"
      ),
      arrow = arrow(
        ends = "last",
        type = "closed",
        angle = 20
      )
    ) + scale_fill_manual(
      name="group",
      values = c(observed = NA, background = "#666666"),
      guide = "none"
    ) +
    scale_colour_manual(
      name="group",
      values = c(background = "#666666", observed = "red"),
      guide = "none"
    ) +
    guides(colour = guide_legend(override.aes = list(
      shape = NA,
      colour = NA,
      fill = c(background =
                 "#666666", observed = "red")
    ))) +
    labs(title = title,
         x="background ratios \n(genes in lemire set / genes not in lemire set)")
  return(list(bg_ors = toplot, obs_or = toplot2, plot = gp))
}

# gather the individual plots
r1 <- enrich_gene_set(obs_genes, lemire_regs, title="Epigenetic regulators")
r2 <- enrich_gene_set(obs_genes, lemire_tfs, title="Transcription factors")
r3 <- enrich_gene_set(obs_genes, lemire_znfs, title="ZNFs")

# save to disc
pdf(fout_plot, width = 16, height=8)
ggarrange(r1$plot, r2$plot, r3$plot,
          ncol=3, common.legend = T, legend="right", labels = "AUTO")
dev.off()

save(file=fout_data, r1, r2, r3)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")

