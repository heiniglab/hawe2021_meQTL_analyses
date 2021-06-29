## annotation of CpG sites with functional genomics data and TFBS predictions

## speed up!
library(parallel)
options(mc.cores=4)

## compute TF affinities with the tRap package
library(tRap)

## load the CpG positions from the bioconductor package
if (!require(FDb.InfiniumMethylation.hg19)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("FDb.InfiniumMethylation.hg19")
}

if (!require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)) {
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
}

## genome sequences
if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Hsapiens.UCSC.hg19")
}


## get 100bp intervals around the CpG sites
cpgs = features(FDb.InfiniumMethylation.hg19)
context = resize(cpgs, 100, fix="center")

compute.affinities <- FALSE
if (compute.affinities) {
  ## get the corresponding sequence
  seq = as.character(getSeq(Hsapiens, context))

  ## get the PWMs
  data(transfac)
  
  matrices = read.transfac("data/current/pwms/znf333.txt")
  
  ## onlly use vertebrate factors

  matrices = matrices[substr(names(matrices), 1, 1) == "V"]

  affinities = mclapply(matrices, function(pwm) affinity(pwm, seq, both.strands=T))

  save(affinities, file="results/current/cpgs_with_affinties.RData")
}

## also annotate CpGs with ChIP-seq binding sites
library(rtracklayer)
library(data.table)


tfbs = import("data/current/tfbs/filPeaks_public.bed")
ann = t(matrix(unlist(strsplit(values(tfbs)[,"name"], ".", fixed=T)), nrow=3))
colnames(ann) = c("geo_id", "TF", "condition")
values(tfbs) = DataFrame(name=values(tfbs)[,"name"], data.frame(ann, stringsAsFactors=F))

## we write out a table with all conditions and select the blood related ones
conditions = t(matrix(unlist(strsplit(unique(values(tfbs)[,"name"]), ".", fixed=T)), nrow=3))
colnames(conditions) = c("geo_id", "TF", "condition")
conditions = conditions[order(conditions[,"condition"]),]
conditions = conditions[,c(1,3)]
conditions = conditions[!duplicated(paste(conditions[,1], conditions[,2])),]

conditions = data.frame(conditions, blood.related=F)

for (term in c("amlpz12_leukemic", "aplpz74_leukemia", "bcell", "bjab", "bl41", "blood", "lcl", "erythroid", "gm", "hbp", "k562", "kasumi", "lymphoblastoid", "mm1s", "p493", "plasma", "sem", "thp1", "u937")) {
  conditions[grep(term, conditions[,2]),"blood.related"] = TRUE
}


selected = tfbs[values(tfbs)[,"condition"] %in% conditions[conditions[,"blood.related"],"condition"]]


## load the encode tfs separately
encode = as.data.frame(fread("data/current/tfbs/wgEncodeRegTfbsClusteredWithCellsV3.bed", header=F))
encode = GRanges(seqnames=encode[,1], ranges=IRanges(encode[,2] + 1, encode[,3]), name=paste("ENCODE", encode[,4], tolower(encode[,6]), sep="."), geo_id="ENCODE", TF=encode[,4], condition=tolower(encode[,6]))

encode.lcl = encode[grep("gm", values(encode)[,"condition"])]
values(encode.lcl)[,"condition"] = "lcl"
encode.k562 = encode[grep("k562", values(encode)[,"condition"])]
values(encode.k562)[,"condition"] = "k562"

selected = c(selected, encode.lcl, encode.k562)
save(selected, file="results/current/tfbs_granges.RData")

chip = paste(values(selected)[,"TF"], values(selected)[,"condition"], sep=".")
chip.exp = unique(chip)



## create an annotation matrix for the CpGs
tfbs.ann = sapply(chip.exp, function(x) overlapsAny(context, selected[chip == x]))
rownames(tfbs.ann) = names(cpgs)

save(tfbs.ann, file="results/current/cpgs_with_chipseq.RData")
save(tfbs.ann, file="results/current/cpgs_with_chipseq_context_100.RData")

chip.exp = data.frame(chip.exp, t(matrix(unlist(strsplit(chip.exp, ".", fixed=T)), nrow=2)))
colnames(chip.exp) = c("experiment", "TF", "conditon")

write.table(chip.exp, file="results/current/chip_experiments.txt", sep="\t", quote=F, row.names=F)


## we also try to merge the information from multiple cell lines for the
## same TFs
ann.by.factor = tapply(as.character(chip.exp$experiment), chip.exp$TF, function(exp.id) apply(tfbs.ann[,exp.id,drop=F], 1, any))
tname = names(ann.by.factor)
ann.by.factor = matrix(unlist(ann.by.factor), ncol=length(ann.by.factor), dimnames=list(names(ann.by.factor[[1]]), tname))
tfbs.ann = ann.by.factor
save(tfbs.ann, file="results/current/cpgs_with_chipseq_context_100_by_tf.RData")


## also perform a variation of this analysis with just the CpG (without context)
## and also larger contexts (500, 1000, 5000, 10000)

sizes = c(2, 10, 500, 1000, 5000, 10000)

for (size in sizes) {
  cat("size:", size, "\n")
  new.context = resize(cpgs, size, fix="center")

  tfbs.ann = sapply(unique(chip), function(x) overlapsAny(new.context, selected[chip == x]))
  rownames(tfbs.ann) = names(cpgs)

  save(tfbs.ann, file=paste0("results/current/cpgs_with_chipseq_context_", size, ".RData"))
}


## annotate with centromere info
cytobands = read.table("data/current/annotation/cytoBand.txt", stringsAsFactors=F)
centromer = cytobands[cytobands[,5] == "acen",]
cstart = tapply(centromer[,2], centromer[,1], min)
cend = tapply(centromer[,3], centromer[,1], max)
centromer = GRanges(seqnames=names(cstart), ranges=IRanges(cstart, cend))

dist = distanceToNearest(cpgs, centromer)

centromere.ann = data.frame(in.centromere=overlapsAny(cpgs, centromer), dist=elementMetadata(dist)[,"distance"])
rownames(centromere.ann) = names(cpgs)


save(centromere.ann, file="results/current/cpgs_with_centromere.RData")
