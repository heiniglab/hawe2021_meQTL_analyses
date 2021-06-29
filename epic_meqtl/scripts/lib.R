# ------------------------------------------------------------------------------
# Library script containing useful methods for the EPIC analysis.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#' Helper to quickly load and annotate ST1 pairs (incl. discovery significance
#' filter of 1e-14)
#'
#' @param fst1 file containing the ST1 pairs
#'
# ------------------------------------------------------------------------------
load_st1_pairs <- function(fst1) {
  read_tsv(fst1) %>% 
  filter(p.eur.discovery < 1e-14) %>%
  mutate(category = case_when(snp.chr == cpg.chr & abs(cpg.pos - snp.pos)>1e6 ~ "longrange",
                              snp.chr == cpg.chr & abs(cpg.pos - snp.pos)<=1e6 ~ "cis",
                              snp.chr != cpg.chr ~ "trans"))
}

#' ------------------------------------------------------------------------------
#' Loads genotype data for the given (snp) ranges.
#'
#' @param fdosage File with dosages to be used (ideally the combined dosages for
#' all SNPs on all chromosomes. Needs to be tabix indexed.
#' @param geno_indiv List of Ids for the genotyped individuals (same order as in
#' the dosage files)
#' @param ranges Ranges in which to look for genotype data. e.g. one range per 
#' SNP
#' @param genotypes_only Whether to return only the genotype information,
#' dropping chr, pos and allele information. Default: T
# ------------------------------------------------------------------------------
load_genotypes <- function(fdosage, geno_indiv, ranges, genotypes_only = T) {
  require(Rsamtools)

  res <- scanTabix(fdosage, param=ranges)

  # scanTabix will also return 'empty' ranges for which no data was found
  res_not_empty <- res[sapply(res, function(x) !length(x)<1)]
  
  if(length(res_not_empty) < 2) {
    collapsed <- paste0(unlist(res_not_empty), "\n")
  } else {
    collapsed <- paste0(unlist(res_not_empty), collapse="\n")
  }
  
  data <- unique(read_tsv(file=collapsed, col_names=FALSE,
                          col_types=cols(X4=col_character(),
                                         X5=col_character())))
 

  if(nrow(data) < 1) {
    warning("No genotypes found for specified ranges")
    return(NULL)
  }

  # TODO: test implementation
  colnames(data) <- c("chr", "snp_id", "pos", "A1", "A2", geno_indiv)
  if(genotypes_only) {
    dat <- t(data[,-c(1:5)])
    colnames(dat) <- data$snp_id
    # indiv x snp matrix incl snp_ids
    return(dat)
  } else {
    # snp x indiv matrix with additional columns (chr, etc.)
    return(data)
  }
}

# ------------------------------------------------------------------------------
# Annotate EPIC meQTL pairs  with the 450k derived sentinels
# ------------------------------------------------------------------------------
annotate_with_sentinels <- function(sentinel_pairs, epic_pairs,
                                    meth_resid, geno_indiv,
                                    r2_cutoff, dist_cutoff,
				    snp_ranges, cpg_ranges,
                                    fdosage, threads) {

  # check r2 for cpgs
  check_cpg_r2 <- function(c1, c2) {
    meth_r2 <- cor(meth_resid[,c1], meth_resid[,c2],
                     use="complete")^2
    if(meth_r2 < r2_cutoff) {
      return(NA)
    } else {
      return(meth_r2)
    }
  }

  # check r2 for snps
  check_snp_r2 <- function(c1, c2) {
    ge <- load_genotypes(fdosage, geno_indiv, snp_ranges[c(c1, c2)])
    if(is.null(ge)) return(NA)
    geno_r2 <- cor(ge[,1], ge[,2], 
                  method="spearman")^2
    if(geno_r2 < r2_cutoff) {
      return(NA)
    } else {
      return(geno_r2)
    }
  }

  sentinel_snps <- (sentinel_pairs %>% pull(snp.sentinel))
  sentinel_cpgs <- (sentinel_pairs %>% pull(cpg.sentinel))

  # already arrange once by p-value (needed later to select
  # top SNP without having to arrange every time...
  epic_pairs <- epic_pairs %>% arrange(`p-value`)

  # we iterate over all CpGs separately
  icpgs <- unique(epic_pairs %>% pull(cpg))

  res <- mclapply((1:length(icpgs)), function(i) {
    ecpg <- icpgs[i]

    epairs_filtered <- epic_pairs %>% filter(cpg == ecpg)

    # list of all associated SNPs
    meqtls <- epairs_filtered %>% pull(SNP)

    # any of the meQTL SNPs a sentinel SNP?
    idx <- which(sentinel_snps %in% meqtls)
    if(length(idx) > 0) {
      # assume sentinel ok, so just check the related CpGs
      ssnp <- sentinel_snps[idx[1]]
      scpgs <- unique(sentinel_cpgs[idx])

      # check distance
      dcpgs <- distance(cpg_ranges[ecpg], cpg_ranges[scpgs])
      didx <- which(!is.na(dcpgs) & dcpgs < dist_cutoff)
      if(length(didx) > 0) {
        for(cg in scpgs[didx]) {
          # check r2
          cpg_r2 <- check_cpg_r2(ecpg, cg)
          if(is.na(cpg_r2)) next;

          # now we have a match
          return(list(snp = ssnp, cpg = cg, 
                      cpg.r2 = cpg_r2, snp.r2 = NA))
        }
      }
    }

    # check cpg dists once
    cpg_dists <- abs(distance(cpg_ranges[ecpg], 
                              cpg_ranges[sentinel_cpgs]))
    cpg_dist_match <- (!is.na(cpg_dists) & cpg_dists < dist_cutoff)

    # for now we check ALL SNPs instead because we had so few matches...
    # select top SNP for cpg (list of SNPs for cpg, already
    # arranged by p-value -> select first
    #esnp <- meqtls[1]

    # check all associated meQTL snps
    for(esnp in meqtls) {

      # get distances
      snp_dists <- abs(distance(snp_ranges[esnp], 
                                snp_ranges[sentinel_snps]))

      # identify sentinel pairs we need to check (matching distances 
      # for both snps and cpgs)
      to_check <- which((!is.na(snp_dists) & snp_dists < dist_cutoff) &
                         cpg_dist_match)

      # we choose a for loop to be able to exit early...
      if(length(to_check) > 0) {
        for(j in to_check) {
          ssnp <- sentinel_snps[j]
          scpg <- sentinel_cpgs[j]
          # debug
          cat(ssnp, scpg, esnp, ecpg, "\n", 
              file="distance_matched_pairs.tsv", sep="\t", append=T)

          # get r2s
          cpg_r2 <- check_cpg_r2(ecpg, scpg)
          if(is.na(cpg_r2)) next;
          snp_r2 <- check_snp_r2(esnp, ssnp)
          if(is.na(snp_r2)) next;

          # debug
          cat(ssnp, scpg, esnp, ecpg, "\n", 
              file="full_matched_pairs.tsv", sep="\t", append=T)

          # we have a match
          return(list(sentinel.snp = ssnp, sentinel.cpg = scpg, 
                      cpg.r2 = cpg_r2, snp.r2 = snp_r2))
        }
      }
    }
  }, mc.cores = threads)
  names(res) <- icpgs
  result <- bind_rows(res, .id = "cpg")

  if(nrow(result) > 0) {
    # match to original data frame
    #left_join(epic_pairs, result, by=c("cpg" , "cpg")) 
    list(epic=epic_pairs, res=result)
  } else {
    warning("No matches found!")
    epic_pairs %>%
      mutate(sentinel.snp = NA_character_, sentinel.cpg = NA_character_, 
             cpg.r2 = NA_real_, snp.r2 = NA_real_)
  }
}

# ------------------------------------------------------------------------------
# Annotate EPIC meQTL pairs  with the 450k derived sentinels
# SNP matches only, returns list of sentinels with matched EPIC CpGs.
# ------------------------------------------------------------------------------
annotate_sentinel_snps_with_epic <- function(sentinel_snps, sentinel_snp_genos, epic_pairs,
                                    geno_indiv, r2_cutoff, dist_cutoff,
                                    snp_ranges, fdosage, threads) {

  # check r2 for snps
  check_snp_r2 <- function(c1, c2, ge) {
    if(!(c1 %in% colnames(ge)) | !(c2 %in% colnames(ge))) {
      return(NA)
    }

    geno_r2 <- cor(ge[,c1], ge[,c2],
                  method="spearman")^2

    if(geno_r2 < r2_cutoff) {
      return(NA)
    } else {
      return(geno_r2)
    }
  }

  # check for each sentinel SNP for matches with the EPIC pairs.
  res <- mclapply(sentinel_snps, function(ssnp) {

    # any of the meQTL SNPs a sentinel SNP?
    direct_matches <- epic_pairs %>% filter(SNP == ssnp)

    # the remainder needs to be checked manually
    pairs_to_check <- anti_join(epic_pairs, direct_matches, 
                                by=c("SNP" = "SNP", "cpg" = "cpg"))

    # filter by distance to sentinel
    dist_matches <- which(abs(distance(snp_ranges[ssnp],
                              snp_ranges[pairs_to_check$SNP])) < dist_cutoff)
    out_cpgs <- c()

    if(length(dist_matches) > 0) {

      # unique distance matched snps
      dist_matched_snps <- pairs_to_check[dist_matches,] %>% pull(SNP) %>% unique

      ge <- load_genotypes(fdosage, geno_indiv, snp_ranges[dist_matched_snps])
      ge <- cbind(sentinel_snp_genos, ge)

      # for the uniquely matched snps, get the CpGs if r2 is OK
      paired_cpgs <- lapply(dist_matched_snps, function(esnp) {

        # get r2
        snp_r2 <- check_snp_r2(esnp, ssnp, ge)
        if(is.na(snp_r2)) return(NULL)

        cpgs <- pairs_to_check %>% filter(SNP == esnp) %>% 
          pull(cpg) %>% unique

        return(cpgs)
      })
      out_cpgs <- unique(unlist(paired_cpgs))
    }
    # if we had direct matches, add them to the collection of cpgs
    if(nrow(direct_matches) > 0) {
      out_cpgs <- unique(c(out_cpgs, direct_matches %>% pull(cpg)))
    }
    paste0(out_cpgs, collapse=",")
  }, mc.cores = threads)

  # return a simple df with sentinel <-> cpgid mapping
  names(res) <- sentinel_snps
  return(res)
}


get_850k_cpg_ranges <- function(subset = NULL) {
  require(tidyverse)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  data(Locations)

  cpg_locs <- Locations %>% as_tibble() %>%
    dplyr::select(chr, start=pos) %>%
    mutate(id=rownames(Locations), end=start+1) %>%
    dplyr::select(id, everything())

  cpg_ranges <- with(cpg_locs, GRanges(chr,IRanges(start,end)))
  names(cpg_ranges) <- cpg_locs$id
  if(!is.null(subset)) {
    return(cpg_ranges[subset])
  } else {
    return(cpg_ranges)
  }
}

# ------------------------------------------------------------------------------
#' recycled method from matthias do determine trans CpGs for a specific sentinel 
#' SNP
# ------------------------------------------------------------------------------
get.trans.cpgs <- function(sentinel, trans.meQTL, cosmo, cpgs=NULL, cosmo.idxs=F) {
  pairs = which(trans.meQTL[,"sentinel.snp"] == sentinel)
  trans.snp.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.snp"], sep=""),
                             ranges=IRanges(trans.meQTL[,"interval.start.snp"],
                                            trans.meQTL[,"interval.end.snp"]))
  trans.cpg.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.cpg"], sep=""),
                             ranges=IRanges(trans.meQTL[,"interval.start.cpg"],
                                            trans.meQTL[,"interval.end.cpg"]))
  pair.snps = GRanges(seqnames=paste("chr", cosmo[,"snp.chr"], sep=""),
                      ranges=IRanges(cosmo[,"snp.pos"], width=1))
  pair.cpgs = GRanges(seqnames=paste("chr", cosmo[,"cpg.chr"], sep=""),
                      ranges=IRanges(cosmo[,"cpg.pos"], width=2))
  pairs = pair.snps %over% trans.snp.ranges[pairs] &
    pair.cpgs %over% trans.cpg.ranges[pairs]
  if(is.null(cpgs) | cosmo.idxs) {
    return(pairs)
  } else {
    is.meQTL = cpgs %over% pair.cpgs[pairs]
    return(is.meQTL)
  }
}


#' -----------------------------------------------------------------------------
#' load TFBS data
#'
#' @param cpgs IDs of CpGs to be annotated with TFBS
#' @param window_size Size of window around CpGs within which to 
#' check for TFBS
#'
#' -----------------------------------------------------------------------------
load_tfbs_annot <- function(cpgs, window_size=100) {

  library(rtracklayer)
  library(data.table)
  tfbs = import("data/tfbs/filPeaks_public.bed")
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

  for (term in c("amlpz12_leukemic", "aplpz74_leukemia", "bcell", "bjab", "bl41", "blood", "lcl", "erythroid",
                 "gm", "hbp", "k562", "kasumi", "lymphoblastoid", "mm1s", "p493", "plasma", "sem", "thp1", "u937")) {
    conditions[grep(term, conditions[,2]),"blood.related"] = TRUE
  }


  selected = tfbs[values(tfbs)[,"condition"] %in% conditions[conditions[,"blood.related"],"condition"]]


  ## load the encode tfs separately
  encode = as.data.frame(fread("data/tfbs/wgEncodeRegTfbsClusteredWithCellsV3.bed", header=F))
  encode = GRanges(seqnames=encode[,1], ranges=IRanges(encode[,2] + 1, encode[,3]),
                   name=paste("ENCODE", encode[,4], tolower(encode[,6]), sep="."),
                   geo_id="ENCODE", TF=encode[,4], condition=tolower(encode[,6]))

  encode.lcl = encode[grep("gm", values(encode)[,"condition"])]
  values(encode.lcl)[,"condition"] = "lcl"
  encode.k562 = encode[grep("k562", values(encode)[,"condition"])]
  values(encode.k562)[,"condition"] = "k562"

  selected = c(selected, encode.lcl, encode.k562)


  chip = paste(values(selected)[,"TF"], values(selected)[,"condition"], sep=".")
  chip.exp = unique(chip)

  cpg_ranges <- get_850k_cpg_ranges()
  context <- resize(cpg_ranges[cpgs], 100, fix="center")
 
  ## create an annotation matrix for the CpGs
  tfbs.ann = sapply(chip.exp, function(x) overlapsAny(context, selected[chip == x]))
  rownames(tfbs.ann) = names(context)
  return(tfbs.ann)
}
