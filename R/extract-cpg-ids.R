# ------------------------------------------------------------------------------
#' Simple script to extract all measured CpG sites from the available KORA data
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------

# get snakemake params
finput_beta <- snakemake@input[["beta"]]
finput_cosmo <- snakemake@input[["cosmo"]]
fout_beta <- snakemake@output[["beta"]]
fout_cosmo_cpg <- snakemake@output[["cosmo_cpg"]]
fout_cosmo_snp <- snakemake@output[["cosmo_snp"]]

# process beta ids
load(finput_beta)
cpgs <- rownames(beta)
cpgs <- cpgs[grepl("^cg", cpgs)]
write.table(file=fout_beta, cpgs, col.names=F, row.names=F, quote=F)

# process cosmo ids
load(finput_cosmo)
cpgs <- unique(as.character(cosmo$cpg))
cpgs <- cpgs[grepl("^cg", cpgs)]
write.table(file=fout_cosmo_cpg, cpgs, col.names=F, row.names=F, quote=F)

snps <- unique(as.character(cosmo$snp))
write.table(file=fout_cosmo_snp, snps, col.names=F, row.names=F, quote=F)
