#!/usr/bin/env Rscript

## look for enrichment of TFBS in clusters of trans meQTLs
source("R/lib.R")

library(parallel)
library(GenomicRanges)
library(optparse)

opt = parse_args(OptionParser(option_list=list(
  make_option("--size", type="integer", default=100),
  make_option("--resample", type="integer", default=0),
  make_option("--cores", type="integer", default=1),
  make_option("--prefix", type="character", default=NULL),
  make_option("--tfbs", type="character", default=NULL),
  make_option("--meQTL", type="character", default=NULL),
  make_option("--array", type="character", default="450k"))))

## set the size of the context region
size = opt$size
nresample = opt$resample
ncores = opt$cores

if (is.null(opt$prefix)) {
  prefix = paste0("results/current/enrichment_chipseq_context_", size, "_resample_", nresample)
} else {
  prefix = opt$prefix
}

## get intervals around the CpG sites
if (opt$array == "450k") {
  require(FDb.InfiniumMethylation.hg19)
  cpgs = features(FDb.InfiniumMethylation.hg19)
}
if (opt$array == "epic") {
  ## problems installing the annotation package.. work around see below
  
  ## BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  ## require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  ## cpgs = features(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  
  ## work around
  library(data.table)
  epic <- fread("data/current/annotation/cpg_annotation_both_arrays.tsv")
  
  epic.cpgs <- GRanges(seqnames=epic[is_epic == TRUE,]$chr, ranges=IRanges(epic[is_epic  == TRUE,]$start, epic[is_epic == TRUE,]$end))
  names(epic.cpgs) <- epic[is_epic == TRUE,]$id
}

context = resize(cpgs, size, fix="center")

## load the annotations
if (is.null(opt$tfbs)) {
  load(paste0("results/current/cpgs_with_chipseq_context_", size, ".RData"))
} else {
  load(paste0(opt$tfbs))
}
  
load("results/current/cpgs_with_centromere.RData")
tfbs.ann = cbind(tfbs.ann, centromere=centromere.ann[,"in.centromere"])

## load the trans clusters
trans.meQTL = read.csv("data/current/meQTLs/transpairs_r02_110117_converted_1MB.txt", sep="\t", stringsAsFactors=F)

if (!"gene.snp" %in% colnames(trans.meQTL)) {
  trans.meQTL = cbind(trans.meQTL, gene.snp=NA)
}

## the pairs will be fetched through the library function to ensure we have the
## exact same pairs!!

if (F) {
  trans.snp.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.snp"], sep=""), ranges=IRanges(trans.meQTL[,"interval.start.snp"], trans.meQTL[,"interval.end.snp"]))

  trans.cpg.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.cpg"], sep=""), ranges=IRanges(trans.meQTL[,"interval.start.cpg"], trans.meQTL[,"interval.end.cpg"]))
}
  
## load the actual pairs to map the intervals back to individual CpGs
load("results/current/trans-cosmopairs_combined_151216.RData")

if (F) {
  pair.snps = GRanges(seqnames=paste("chr", cosmo[,"snp.chr"], sep=""), ranges=IRanges(cosmo[,"snp.pos"], width=1))
  pair.cpgs = GRanges(seqnames=paste("chr", cosmo[,"cpg.chr"], sep=""), ranges=IRanges(cosmo[,"cpg.pos"], width=2))
}

## also perform a test based on resampling of cpgs from a matched background
## load the msd object that has mean and sd of the beta values for all cpgs
load("data/current/meQTLs/msd_lolipop.RData")

## make sure that we have the same order as in the coordinates
msd = msd[match(names(cpgs), rownames(msd)),]
rownames(msd) = names(cpgs)

res = NULL
failed = NULL
for (sentinel in unique(trans.meQTL$sentinel.snp)) {
  ## sentinel = "rs60626639" # CTCF
  
  cat(sentinel, "\n")

  pairs = which(trans.meQTL[,"sentinel.snp"] == sentinel)

  if (length(pairs) < 5) {
    next
  }

  if (F) {
    pairs = pair.snps %over% trans.snp.ranges[pairs] & pair.cpgs %over% trans.cpg.ranges[pairs]
  }
  
  is.meQTL = factor(get.trans.cpgs(sentinel, trans.meQTL, cosmo, cpgs), levels=c(FALSE, TRUE))

  if (nresample > 0) {
    ## get a matching background cpgs for each cpg with meQTL
    matched = lapply(names(cpgs)[is.meQTL == TRUE], get.matched.cpgs, msd)
  }
  
  ## this is the regular loop
  ## for (tf in colnames(tfbs.ann)) {

  ## we would like to speed up things
  rfun <- function(tf) {
    cat(".")
    is.bound = factor(tfbs.ann[,tf], levels=c(FALSE, TRUE))
    
    tab = table(is.meQTL=is.meQTL, is.bound=is.bound)
    test = fisher.test(tab)
    expected = chisq.test(tab)$expected

    set.seed(match(tf, colnames(tfbs.ann)))

    if (nresample > 0) {
      alternative = c("less", "greater")[as.numeric(test$estimate > 1) + 1]
      emp = empirical.enrichment(is.meQTL, is.bound, matched, alternative=alternative, n.resample=nresample, over=0.05)
    } else {
      emp = NA
    }
    
    tab = as.data.frame(tab)
    n = paste(colnames(tab)[1], as.character(tab[,1]), colnames(tab)[2], as.character(tab[,2]), sep=":")
    
    this = data.frame(sentinel, tf, test$estimate, p=test$p.value, empirical.p=emp, t(tab[,"Freq"]), t(as.numeric(expected)))
    colnames(this)[6:ncol(this)] = c(paste("obs", n, sep="."), paste("exp", n, sep="."))
    return(this)
  }
  ## missing = 1:ncol(tfbs.ann)
  ## ans = list()
  ## for (i in 1:10) {
  ##   ans[missing] = mclapply(colnames(tfbs.ann)[missing], rfun, mc.cores=10)
  ##   missing = which(!sapply(ans, inherits, "data.frame")
  ##   if (length(missing) == 0)
  ##     break
  ##   cat("missing:", length(missing), "\n")
  ## }

  ans = mclapply(colnames(tfbs.ann), rfun, mc.cores=ncores)
  
  for (i in 1:length(ans)) {
    if (inherits(ans[[i]], "data.frame")) {
      res = rbind(res, ans[[i]])
    } else {
      failed = rbind(failed, data.frame(sentinel, tf=colnames(tfbs.ann)[i]))
    }
  }
  cat("\n")
}
  
res = cbind(res, q=p.adjust(res$p, "BH"), empirical.q=p.adjust(res$empirical.p, "BH"))

## for a conservative analysis also take the max of empirical and theoretical P-value
res = cbind(res, max.q=pmax(res$q, res$empirical.q))

write.table(res, paste0(prefix, ".txt"), sep="\t", quote=F, row.names=F)

## make the results a bit more accessible
res = with(res, cbind(res,
  tf.symbol=sapply(strsplit(as.character(tf), ".", fixed=T), "[", 1),
  condition=sapply(strsplit(as.character(tf), ".", fixed=T), "[", 2),
  nCpG.meQTL=obs.is.meQTL.TRUE.is.bound.FALSE + obs.is.meQTL.TRUE.is.bound.TRUE,
  nCpG.no.meQTL = obs.is.meQTL.FALSE.is.bound.FALSE +
    obs.is.meQTL.FALSE.is.bound.TRUE))

res = with(res, cbind(res,
  pct.bound.meQTL=obs.is.meQTL.TRUE.is.bound.TRUE / nCpG.meQTL,
  pct.bound.no.meQTL=obs.is.meQTL.FALSE.is.bound.TRUE / nCpG.no.meQTL))

  
write.table(res, paste0(prefix, ".txt"), sep="\t", quote=F, row.names=F)

sig = res[res$q < 0.05,]

trans.ann = trans.meQTL[!duplicated(trans.meQTL$sentinel.snp), c("sentinel.snp", "gene.snp")]

merged = merge(trans.ann, sig, by.x="sentinel.snp", by.y="sentinel", suffixes=c(".qtl", ".tfbs"))
write.table(merged, file=paste0(prefix, "_sig.txt"), sep="\t", quote=F, row.names=F)

## moved heatmap code to lib.R


pdf(file=paste0(prefix, "_heatmap.pdf"), width=8.27, height=15)
print(tfbs.heatmap(res, "q", 0.05, prune=FALSE, max.or=100, min.or=0.01))
dev.off()

pdf(file=paste0(prefix, "_heatmap_pruned.pdf"), width=8.27, height=15)
print(tfbs.heatmap(res, "q", 0.05, prune=TRUE, max.or=100, min.or=0.01))
dev.off()

pdf(file=paste0(prefix, "_empirical_heatmap.pdf"), width=8.27, height=15)
print(tfbs.heatmap(res, "empirical.q", 0.05, prune=FALSE, max.or=100, min.or=0.01))
dev.off()

pdf(file=paste0(prefix, "_empirical_heatmap_pruned.pdf"), width=8.27, height=15)
print(tfbs.heatmap(res, "empirical.q", 0.05, prune=TRUE, max.or=100, min.or=0.01))
dev.off()

pdf(file=paste0(prefix, "_maxq_heatmap.pdf"), width=8.27, height=15)
print(tfbs.heatmap(res, "max.q", 0.05, prune=FALSE, max.or=100, min.or=0.01))
dev.off()

pdf(file=paste0(prefix, "_maxq_heatmap_pruned.pdf"), width=8.27, height=15)
print(tfbs.heatmap(res, "max.q", 0.05, prune=TRUE, max.or=100, min.or=0.01))
dev.off()


if (F) {
## plot as a network graph
library(graph)
library(Rgraphviz)

nodes = union(as.character(sig$sentinel), as.character(sig$tf))
g = graphNEL(nodes=nodes)
g = addEdge(as.character(sig$sentinel), as.character(sig$tf), g)




## also analyse the distace distribution of Cpgs to the centromere
centro = NULL
for (sentinel in unique(trans.meQTL$sentinel.snp)) {
  cat(sentinel, "\n")

  pairs = which(trans.meQTL[,"sentinel.snp"] == sentinel)

  if (length(pairs) < 5) {
    next
  }
  
  pairs = pair.snps %over% trans.snp.ranges[pairs] & pair.cpgs %over% trans.cpg.ranges[pairs]
  
  is.meQTL = factor(cpgs %over% pair.cpgs[pairs], levels=c(FALSE, TRUE))

  ## get a matching background cpgs for each cpg with meQTL
  ## matched = lapply(names(cpgs)[is.meQTL == TRUE], get.matched.cpgs, msd)
  
  w1 = wilcox.test(centromere.ann[,"dist"] ~ is.meQTL)
  means = tapply(centromere.ann[,"dist"], is.meQTL, mean)
  centro = rbind(centro, data.frame(sentinel, estimate=w1$statistic, p=w1$p.value, mean.dist.meQTL=means["TRUE"], mean.dist.no.meQTL=means["FALSE"]))
}

}

if (F) {
## specific CTCF analysis

## in the first round we observed enrichment for CTCF and members of the
## cohesin complex so we define new binding site sets:
## cohesin:
cohesin = with(as.data.frame(tfbs.ann),
  (RAD21.lcl | RAD21.k562) &
  (SMC3.lcl | SMC3.k562) &
  (SMC1A.lcl_hls554p | SMC1A.lcl_bcbl1))
CTCF = with(as.data.frame(tfbs.ann),
  CTCF.lcl | CTCF.k562 | CTCF.bl41 | CTCF.lcl_bcbl1)
CTCF.1.cohesin.1 = cohesin & CTCF
CTCF.1.cohesin.0 = !cohesin & CTCF
CTCF.0.cohesin.1 = cohesin & !CTCF

tfbs.ann = cbind(tfbs.ann, cohesin, CTCF, CTCF.1.cohesin.1, CTCF.1.cohesin.0, CTCF.0.cohesin.1)

## we do a similar analysis for the other potential cofactors
CTCFL = tfbs.ann[,"CTCFL.k562"]
CTCF.1.CTCFL.1 = CTCFL & CTCF
CTCF.1.CTCFL.0 = !CTCFL & CTCF
CTCF.0.CTCFL.1 = CTCFL & !CTCF

tfbs.ann = cbind(tfbs.ann, CTCFL, CTCF.1.CTCFL.1, CTCF.1.CTCFL.0, CTCF.0.CTCFL.1)

ZNF143 = with(as.data.frame(tfbs.ann), ZNF143.k562 | ZNF143.lcl)
CTCF.1.ZNF143.1 = ZNF143 & CTCF
CTCF.1.ZNF143.0 = !ZNF143 & CTCF
CTCF.0.ZNF143.1 = ZNF143 & !CTCF

tfbs.ann = cbind(tfbs.ann, ZNF143, CTCF.1.ZNF143.1, CTCF.1.ZNF143.0, CTCF.0.ZNF143.1)

ARID3A = tfbs.ann[,"ARID3A.k562"]
CTCF.1.ARID3A.1 = ARID3A & CTCF
CTCF.1.ARID3A.0 = !ARID3A & CTCF
CTCF.0.ARID3A.1 = ARID3A & !CTCF

tfbs.ann = cbind(tfbs.ann, ARID3A, CTCF.1.ARID3A.1, CTCF.1.ARID3A.0, CTCF.0.ARID3A.1)

MAFK = tfbs.ann[,"MAFK.k562"]
CTCF.1.MAFK.1 = MAFK & CTCF
CTCF.1.MAFK.0 = !MAFK & CTCF
CTCF.0.MAFK.1 = MAFK & !CTCF

tfbs.ann = cbind(tfbs.ann, MAFK, CTCF.1.MAFK.1, CTCF.1.MAFK.0, CTCF.0.MAFK.1)

MAZ = tfbs.ann[,"MAZ.k562"]
CTCF.1.MAZ.1 = MAZ & CTCF
CTCF.1.MAZ.0 = !MAZ & CTCF
CTCF.0.MAZ.1 = MAZ & !CTCF

tfbs.ann = cbind(tfbs.ann, MAZ, CTCF.1.MAZ.1, CTCF.1.MAZ.0, CTCF.0.MAZ.1)

RCOR1 = tfbs.ann[,"RCOR1.k562"]
CTCF.1.RCOR1.1 = RCOR1 & CTCF
CTCF.1.RCOR1.0 = !RCOR1 & CTCF
CTCF.0.RCOR1.1 = RCOR1 & !CTCF

tfbs.ann = cbind(tfbs.ann, RCOR1, CTCF.1.RCOR1.1, CTCF.1.RCOR1.0, CTCF.0.RCOR1.1)

cofactor = CTCFL | ZNF143 | ARID3A | MAFK | MAZ | RCOR1
CTCF.1.cofactor.1 = cofactor & CTCF
CTCF.1.cofactor.0 = !cofactor & CTCF
CTCF.0.cofactor.1 = cofactor & !CTCF

tfbs.ann = cbind(tfbs.ann, cofactor, CTCF.1.cofactor.1, CTCF.1.cofactor.0, CTCF.0.cofactor.1)

sentinel = trans.meQTL$sentinel.snp[trans.meQTL$gene.snp == "CTCF"][1]

pairs = which(trans.meQTL[,"sentinel.snp"] == sentinel)
pairs = pair.snps %over% trans.snp.ranges[pairs] & pair.cpgs %over% trans.cpg.ranges[pairs]
  
is.meQTL = factor(cpgs %over% pair.cpgs[pairs], levels=c(FALSE, TRUE))

## get a matching background cpgs for each cpg with meQTL
matched = lapply(names(cpgs)[is.meQTL == TRUE], get.matched.cpgs, msd)

bg = unlist(lapply(matched, sample, 10))

save(list=c("tfbs.ann", "sentinel", "is.meQTL", "matched"), file="ctcf.RData")

selected = c(which(is.meQTL == TRUE), bg)
is.meQTL = is.meQTL[selected]


tfs = c("cohesin", "RAD21.lcl", "SMC3.lcl", "SMC1A.lcl_hls554p", "CTCF", "CTCF.1.cohesin.1", "CTCF.1.cohesin.0", "CTCF.0.cohesin.1", "CTCFL", "CTCF.1.CTCFL.1", "CTCF.1.CTCFL.0", "CTCF.0.CTCFL.1", "ZNF143", "CTCF.1.ZNF143.1", "CTCF.1.ZNF143.0", "CTCF.0.ZNF143.1", "ARID3A", "CTCF.1.ARID3A.1", "CTCF.1.ARID3A.0", "CTCF.0.ARID3A.1", "MAFK", "CTCF.1.MAFK.1", "CTCF.1.MAFK.0", "CTCF.0.MAFK.1", "MAZ", "CTCF.1.MAZ.1", "CTCF.1.MAZ.0", "CTCF.0.MAZ.1", "RCOR1", "CTCF.1.RCOR1.1", "CTCF.1.RCOR1.0", "CTCF.0.RCOR1.1", "cofactor", "CTCF.1.cofactor.1", "CTCF.1.cofactor.0", "CTCF.0.cofactor.1")


ctcf.res = NULL
for (tf in tfs) {
  cat(".")
  is.bound = factor(tfbs.ann[selected,tf], levels=c(FALSE, TRUE))
  
  tab = table(is.meQTL=is.meQTL, is.bound=is.bound)
  test = fisher.test(tab)
  expected = chisq.test(tab)$expected
  
  ## set.seed(match(tf, colnames(tfbs.ann)))
  ## emp = empirical.enrichment(is.meQTL, is.bound, matched)
  emp = c(p.lt=NA, p.gt=NA)
  names(emp) = paste("emp", names(emp), sep=".")
  
  tab = as.data.frame(tab)
  n = paste(colnames(tab)[1], as.character(tab[,1]), colnames(tab)[2], as.character(tab[,2]), sep=":")
  
  this = data.frame(sentinel, tf, test$estimate, p=test$p.value, t(emp), t(tab[,"Freq"]), t(as.numeric(expected)))
  colnames(this)[7:ncol(this)] = c(paste("obs", n), paste("exp", n))
  ctcf.res = rbind(ctcf.res, this)
  cat("\n")
}

ctcf.res = cbind(ctcf.res, 
  pct.bound.meQTL = ctcf.res[,"obs is.meQTL:TRUE:is.bound:TRUE"] /
  (ctcf.res[,"obs is.meQTL:TRUE:is.bound:FALSE"] +
   ctcf.res[,"obs is.meQTL:TRUE:is.bound:TRUE"]),
  pct.bound.no.meQTL = ctcf.res[,"obs is.meQTL:FALSE:is.bound:TRUE"] /
  (ctcf.res[,"obs is.meQTL:FALSE:is.bound:FALSE"] +
   ctcf.res[,"obs is.meQTL:FALSE:is.bound:TRUE"]))

write.table(ctcf.res, file=paste0("results/current/ctcf-cofactors-matched-context-", size, ".txt"), sep="\t", quote=F, row.names=F)


library(ggplot2)
library(reshape)
plot.data = melt(ctcf.res[,c("tf", "pct.bound.meQTL", "pct.bound.no.meQTL")])
colnames(plot.data) = c("tf", "meQTL", "pct.bound")
plot.data$meQTL = factor(plot.data$meQTL, levels=c("pct.bound.meQTL", "pct.bound.no.meQTL"), labels=c("meQTL", "no meQTL"))

pdf(file=paste0("results/current/ctcf-cofactors-matched-context-", size, ".pdf"), width=14)
ggplot(aes(tf, pct.bound), data=plot.data) + xlab("Transcription factor") + ylab("Percent CpG bound by TF") + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_grid(. ~ meQTL)
dev.off()

}
