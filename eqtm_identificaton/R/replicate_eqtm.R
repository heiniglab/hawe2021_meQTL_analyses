#' -----------------------------------------------------------------------------
#' Check replication of eQTM results from KORA and LOLIPOP SA
#' Script assumes input files sorted by CpG and Gene names
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Jun  12 07:47:37 2020
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

# total number of lines (w/o header) in files is: 5,738,181,750

# location files
fcpg_pos <- snakemake@input$cpg_positions
fgene_pos <- snakemake@input$gene_positions

# load position annotations immediately
load(fcpg_pos) # loads 'cpg.pos'
cpg_pos <- cpg.pos
gene_pos <- read_delim(fgene_pos, delim = ",") %>%
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

# kora
feur <- snakemake@input$eur
# lolipop
fsa <- snakemake@input$sa

# output file for meta associations
fout_kora <- snakemake@output$kora
fout_lolipop <- snakemake@output$lolipop
fout_counts <- snakemake@output$counts

# total number of pairs
total_pairs <- 5738181750

# FDR cutoff
fdr_cut <- 0.05
# bonferroni cutoff
bonf_cut <- 0.05 / total_pairs
# replication cutoff
repl_cut <- 0.05

threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Processing.")
# ------------------------------------------------------------------------------

# alternative method using a filehandle and convertig to data.table
read_next <- function(fhandle, size) {
  # this seems to be the best option for now.
  # fread() throws an error (pot bug), and read_tsv() does not allow
  # to seek to the end of files (also reports error)
  tmp <- readLines(fhandle, n = size)
  fread(paste(tmp, collapse="\n"), 
        col.names = c("cpg", "gene", "beta", "tstat", "pvalue", "P.BH"),
        nThread = threads)
}

feur_handle <- file(feur, open="r")
fsa_handle <- file(fsa, open="r")

# read one line to avoid the header
h <- readLines(feur_handle,1)
h <- readLines(fsa_handle,1)

# prepare main loop;
read <- T
size <- 1e7
i <- 1

ty <- Sys.time()

# prepare the replication counters
kora_sign_fdr <- 0
lol_sign_fdr <- 0
kora_sign_bonf <- 0
lol_sign_bonf <- 0

kora_rep_in_lol_fdr <- 0
lol_rep_in_kora_fdr <- 0
kora_rep_in_lol_bonf <- 0
lol_rep_in_kora_bonf <- 0

# we need individual counters for the specific categories
counters <- c(kora_sign_fdr=kora_sign_fdr,
              lol_sign_fdr=lol_sign_fdr,
              kora_sign_bonf=kora_sign_bonf,
              lol_sign_bonf=lol_sign_bonf,
              kora_rep_in_lol_fdr=kora_rep_in_lol_fdr,
              lol_rep_in_kora_fdr=lol_rep_in_kora_fdr,
              kora_rep_in_lol_bonf=kora_rep_in_lol_bonf,
              lol_rep_in_kora_bonf=lol_rep_in_kora_bonf)

all_counters <- list(cis = counters,
                     longrange = counters,
                     trans = counters)

update_counters <- function(all_counters, eur, sa) {
  
  categories <- names(all_counters)
  
  all_counters <- lapply(categories, function(ca) {
    
    cnts <- all_counters[[ca]]
    
    eur_sub <- filter(eur, category == ca)
    sa_sub <- filter(sa, category == ca)
    
    # get info whether beta-signs match
    beta_match <- (sign(eur_sub$beta) == sign(sa_sub$beta))
  
    # get the counts
    cnts["kora_sign_fdr"] <- cnts["kora_sign_fdr"] + sum(eur_sub$sign_fdr)
    cnts["lol_sign_fdr"] <- cnts["lol_sign_fdr"] + sum(sa_sub$sign_fdr)
    cnts["kora_sign_bonf"] <- cnts["kora_sign_bonf"] + sum(eur_sub$sign_bonf)
    cnts["lol_sign_bonf"] <- cnts["lol_sign_bonf"] + sum(sa_sub$sign_bonf)
    
    # check replication, also check direction of effect
    cnts["kora_rep_in_lol_fdr"] <- cnts["kora_rep_in_lol_fdr"] + 
      sum(eur_sub$sign_fdr & sa_sub$sign_rep & beta_match)
    cnts["lol_rep_in_kora_fdr"] <- cnts["lol_rep_in_kora_fdr"] + 
      sum(sa_sub$sign_fdr & eur_sub$sign_rep & beta_match)
    cnts["kora_rep_in_lol_bonf"] <- cnts["kora_rep_in_lol_bonf"] + 
      sum(eur_sub$sign_bonf & sa_sub$sign_rep & beta_match)
    cnts["lol_rep_in_kora_bonf"] <- cnts["lol_rep_in_kora_bonf"] + 
      sum(sa_sub$sign_bonf & eur_sub$sign_rep & beta_match)
    
    cnts
  })
  names(all_counters) <- categories
  return(all_counters)
}

print("Starting loop...")
while(read) {
  if(i %% 10 == 0) {
    print(paste0("Progress: ", round(i*size / total_pairs * 100, digits=2)))
    print("Time elapsed:")
    print(Sys.time() - ty)
    gc(full=T)
  }
  # get data chunks
  print("Loading data...")
  eur <- read_next(feur_handle, size) %>%
    mutate(sign_bonf = pvalue < bonf_cut,
           sign_fdr = P.BH < fdr_cut,
           sign_rep = pvalue < repl_cut)
  sa <- read_next(fsa_handle, size) %>%
    mutate(sign_bonf = pvalue < bonf_cut,
           sign_fdr = P.BH < fdr_cut,
           sign_rep = pvalue < repl_cut)
  
  # get distances
  dists <- distance(cpg_ranges[eur$cpg], gene_ranges[eur$gene])
  
  eur <- mutate(eur,
                cpg_chr = as.character(seqnames(cpg_ranges[cpg])),
                cpg_pos = start(cpg_ranges[cpg]),
                gene_chr = as.character(seqnames(gene_ranges[gene])),
                gene_start = start(gene_ranges[gene]),
                gene_end = end(gene_ranges[gene]),
                category = case_when(is.na(dists) ~ "trans",
                               !is.na(dists) & abs(dists) > 1e6 ~ "longrange",
                               !is.na(dists) & abs(dists) <= 1e6 ~ "cis"))
  sa <- mutate(sa,
               cpg_chr = as.character(seqnames(cpg_ranges[cpg])),
               cpg_pos = start(cpg_ranges[cpg]),
               gene_chr = as.character(seqnames(gene_ranges[gene])),
               gene_start = start(gene_ranges[gene]),
               gene_end = end(gene_ranges[gene]),
               category = case_when(is.na(dists) ~ "trans",
                                    !is.na(dists) & abs(dists) > 1e6 ~ "longrange",
                                    !is.na(dists) & abs(dists) <= 1e6 ~ "cis"))
  
  print("Updating counts..")
  all_counters <- update_counters(all_counters, eur, sa)
  
  print("Saving FDR results.")
  if(i == 1) {
    fwrite(eur %>% filter(sign_fdr) %>% 
             dplyr::select(-sign_fdr, -sign_bonf, -sign_rep), 
           fout_kora, append = F, sep="\t",
           nThread = threads)
    fwrite(sa %>% filter(sign_fdr) %>% 
             dplyr::select(-sign_fdr, -sign_bonf, -sign_rep), 
           fout_lolipop, append = F, sep="\t",
           nThread = threads)
  } else {
    fwrite(eur %>% filter(sign_fdr) %>% 
             dplyr::select(-sign_fdr, -sign_bonf, -sign_rep), 
           fout_kora, append = T, sep="\t",
           nThread = threads)
    fwrite(sa %>% filter(sign_fdr) %>% 
             dplyr::select(-sign_fdr, -sign_bonf, -sign_rep), 
           fout_lolipop, append = T, sep="\t",
           nThread = threads)
  }
  
  i <- i + 1
  
  if(nrow(eur) < size) read <- F
}
print("Loop done.")
close(feur_handle)
close(fsa_handle)

print("Total time elapsed:")
print(Sys.time() - ty)

saveRDS(all_counters, fout_counts)

tmp <- lapply(c("cis", "longrange", "trans"), function(ca) {
  
  cnts <- as.list(all_counters[[ca]])
  
  print(paste0(ca, " --------------------------------------------------------"))
  # report numbers
  print("Significant in KORA (bonf/FDR):")
  print(paste0(cnts$kora_sign_bonf, " / ", cnts$kora_sign_fdr))
  print("Replication in LOLIPOP (bonf/FDR):")
  print(paste0(cnts$kora_rep_in_lol_bonf, " / ", cnts$kora_rep_in_lol_fdr))
  
  print("Significant in LOLIPOP (bonf/FDR):")
  print(paste0(cnts$lol_sign_bonf, " / ", cnts$lol_sign_fdr))
  print("Replication in KORA (bonf/FDR):")
  print(paste0(cnts$lol_rep_in_kora_bonf, " / ", cnts$lol_rep_in_kora_fdr))
  
})

# ------------------------------------------------------------------------------
print("Done.SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
