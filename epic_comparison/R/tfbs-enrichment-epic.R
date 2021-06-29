#' -----------------------------------------------------------------------------
#' Perform the TFBS enrichment analysis for the EPIC discovery meQTLs and 
#' compare with 450k based enrichment.
#'
#' @author Matthias Heinig
#'
#' @date Mon Mar  9 15:03:59 2020
#' -----------------------------------------------------------------------------

# snakemake logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
source("R/lib.R")

library(parallel)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(parallel)

# ------------------------------------------------------------------------------
# Prepare params

# params
ncores <- snakemake@threads
nresample <- as.numeric(snakemake@wildcards$nresample)
cpg_subset <- snakemake@wildcards$cpg_subset

print("Number of resamplings:")
print(nresample)

# output files
fout_results <- snakemake@output$result
fout_tfbs_annot <- snakemake@output$tfbs_annot

# input files
fcpg_annot <- snakemake@input$cpg_annot
fmeqtls <- snakemake@input$meqtls
ftfbs_granges <- snakemake@input$tfbs_granges

# ------------------------------------------------------------------------------
# Load and prep data

## load epic data
epic <- fread(fcpg_annot)
# only subset where we also want a different background
if("450kwith450kBG" %in% cpg_subset) {
  epic <- subset(epic, is_450k & !is.na(chr))
} else {
  epic <- subset(epic, is_epic)  
}

cpgs <- with(epic, GRanges(chr, 
                           IRanges(start,
                                   end)))
names(cpgs) <- epic$id

context = resize(cpgs, 100, fix="center")

# ------------------------------------------------------------------------------
## load tfbs annotations

load(file=ftfbs_granges) ## load "selected"

chip = paste(values(selected)[,"TF"], values(selected)[,"condition"], sep=".")
chip.exp = unique(chip)

tfbs.ann = sapply(chip.exp, function(x) overlapsAny(context, selected[chip == x]))
rownames(tfbs.ann) = names(cpgs)

save(tfbs.ann, file=fout_tfbs_annot)

# ------------------------------------------------------------------------------
# Load meQTLs. here we only have sentinel SNPs in the table
meQTLs <- read.csv(fmeqtls,
                   sep="\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# define which CpGs to use: either 450k CpGs or 450k + independent EPIC

if("450konly" %in% cpg_subset | "450kwith450kBG" %in% cpg_subset) {
  meQTLs <- subset(meQTLs, on_450k)
} else if("450kwithEPIC" %in% cpg_subset) {
  meQTLs <- subset(meQTLs, on_450k | !has_450k_proxy)
} else if("EPICFull" %in% cpg_subset) {
  meQTLs <- meQTLs
} else if("EPICnoNovel" %in% cpg_subset) {
  meQTLs <- subset(meQTLs, on_450k | has_450k_proxy)
}else {
  stop(paste0("CpG subset not supported: ", cpg_subset))
}

# ------------------------------------------------------------------------------
# get the msd data.frame for all cpgs on the EPIC array

msd <- as.data.frame(epic)
rownames(msd) <- msd$id

# the only subset where we need a different (450k) background
if("450kwith450kBG" %in% cpg_subset) {
  msd <- msd[,c("mean.450k", "sd.450k")]
} else {
  msd <- msd[,c("mean.epic", "sd.epic")]
}

colnames(msd) <- c("meanBeta", "sdBeta")

# ------------------------------------------------------------------------------
# Perform analysis

res = NULL
failed = NULL
for (sentinel in unique(meQTLs$SNP)) {
  ## sentinel = "rs60626639" # CTCF
  
  cat(sentinel, "\n")
  
  pairs = which(meQTLs[,"SNP"] == sentinel)
  
  if (length(pairs) < 5) {
    next
  }
  
  is.meQTL = factor(names(cpgs) %in% meQTLs[pairs,"cpg"] , levels=c(FALSE, TRUE))
  
  if (nresample > 0) {
    ## get a matching background cpgs for each cpg with meQTL
    matched = mclapply(names(cpgs)[is.meQTL == TRUE], get.matched.cpgs, msd, 
                       mc.cores = ncores)
  }
  
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
  ans = mclapply(colnames(tfbs.ann), rfun, 
                 mc.cores=ncores)
  
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

## make the results a bit more accessible
res = with(res, cbind(res,
                      tf.symbol=sapply(strsplit(as.character(tf), ".", fixed=T), "[", 1),
                      condition=sapply(strsplit(as.character(tf), ".", fixed=T), "[", 2),
                      nCpG.meQTL=`obs.is.meQTL:TRUE:is.bound:FALSE` + `obs.is.meQTL:TRUE:is.bound:TRUE`,
                      nCpG.no.meQTL = `obs.is.meQTL:FALSE:is.bound:FALSE` +
                        `obs.is.meQTL:FALSE:is.bound:TRUE`))

res = with(res, cbind(res,
                      pct.bound.meQTL=`obs.is.meQTL:TRUE:is.bound:TRUE` / nCpG.meQTL,
                      pct.bound.no.meQTL=`obs.is.meQTL:FALSE:is.bound:TRUE` / nCpG.no.meQTL))


# all done, write output
write.table(res, fout_results, sep="\t", quote=F, row.names=F)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

