#' -----------------------------------------------------------------------------
#' Methods for the snp-cpg annotation overlap enrichment analysis
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jun 11 11:20:21 2019
#' -----------------------------------------------------------------------------

#' -----------------------------------------------------------------------------
#' Check the overlap in states for a set of CpG-SNP pairs (single CpG)
#'
#' @param cpg The cpg id
#' @param snps The SNP ids
#' @param cpg_states List of chromHMM states for CpGs (named by CpGs)
#' @param snp_states List of chromHMM states for SNPs (named by SNPs)
#' @param enhancer Vector of chromHMM states considered as enhancer states
#' @param promoter Vector of chromHMM states considered as promoter states
#'
#' -----------------------------------------------------------------------------
get_overlap_region <- function(cpg, snps, cpg_states, snp_states,
                               enhancer, promoter) {
  res <- c()
  cpg_state <- cpg_states[cpg]
  snp_states_sub <- snp_states[snps]

  enhancer_snp <- names(snp_states_sub[snp_states_sub %in% enhancer])
  promoter_snp <- names(snp_states_sub[snp_states_sub %in% promoter])

  if(is.na(cpg_state) |
     all(is.na(snp_states))) {
    res <- rep(NA, 4)
  } else {
    if (cpg_state %in% enhancer &
        length(enhancer_snp) > 0) {
      res <- 1
    } else {
      res <- 0
    }
    if (cpg_state %in% enhancer &
        length(promoter_snp) > 0) {
      res <- c(res, 1)
    } else {
      res <- c(res, 0)
    }
    if (cpg_state %in% promoter &
        length(enhancer_snp) > 0) {
      res <- c(res, 1)
    } else {
      res <- c(res, 0)
    }
    if (cpg_state %in% promoter &
        length(promoter_snp) > 0) {
      res <- c(res, 1)
    } else {
      res <- c(res, 0)
    }
  }

  names(res) <- c("enhancer", "enh_prom",
                  "prom_enh", "promoter")

  # get additional information on SNP (remember one example or set NA)
  if(length(enhancer_snp) > 0) {
    enhancer_snp <- enhancer_snp[1]
  } else {
    enhancer_snp <- NA_character_
  }
  if(length(enhancer_snp) > 0) {
    promoter_snp <- promoter_snp[1]
  } else {
    promoter_snp <- NA_character_
  }

  result <- list(overlaps = res,
                 enhancer_snp = enhancer_snp,
                 promoter_snp = promoter_snp)
  return(result)
}

#' -----------------------------------------------------------------------------
#' Gets a matched background region for the provided CpG and its associated SNPs
#' -----------------------------------------------------------------------------
get_matched_region <- function(cpg,
                               snps,
                               cosmo_cis,
                               snp_freq_by_chr,
                               cpg_msd,
                               snps_background,
                               cpgs_background,
                               threads = 1) {

  require(parallel)

  # allow +-1000bp distance mismatch
  dist_eps <- 1000

  # allow 0.05 MAF difference
  maf_eps <- 0.05

  # distances between CpG and associated SNPs
  dists <- abs(cosmo_cis[match(snps, cosmo_cis$snp), "snp.pos"] -
                 cosmo_cis[match(cpg, cosmo_cis$cpg), "cpg.pos"])
  names(dists) <- snps

  # get the mafs of the QTL SNPs
  qtl_snp_chr <- unique(cosmo_cis[match(snps, cosmo_cis$snp), "snp.chr"])
  mafs <- snp_freq_by_chr[[qtl_snp_chr]][match(snps,
                                               snp_freq_by_chr[[qtl_snp_chr]]$rsid), "maf"]
  names(mafs) <- snps

  matched_cpgs <- rownames(cpg_msd[get.matched.cpgs(cpg, cpg_msd), ])
  matched_cpgs <-
    sample(intersect(cpgs_background$cpg, matched_cpgs)) # sample cpgs

  # get entity position information to match by distance
  cpg_pos <- cpgs_background[matched_cpgs, , drop = F]

  # no matched CpGs
  if(nrow(cpg_pos) < 1) {
    return(NA)
  }

  best_match <- NA
  best_match_nas <- length(dists)

  print(paste0("Checking ", nrow(cpg_pos), " cpgs."))

  # check all matched CpGs
  for (i in 1:nrow(cpg_pos)) {
    #print(i)
    icpg <- cpg_pos[i, "cpg"]
    ichr <- cpg_pos[i, "chr"]
    ipos <- cpg_pos[i, "pos"]
  #  print(ichr)
    # get all SNPs, which are within the max distance
    snp_freq_sub <- subset(snp_freq_by_chr[[ichr]],
                             abs(ipos - pos) <= max(dists) + dist_eps &
                             rsid %in% snps_background$snp)
    snp_freq_sub$dist <- abs(snp_freq_sub$pos - ipos)

    # only proceed if we have as least as many snps as we need within the window
    if(nrow(snp_freq_sub) < length(dists)) next

    # for each SNP (named distances) get a matched background SNP such that each
    # distance is reflected once
    matched_snps <- mclapply(1:length(dists), function(j) {

      d <- dists[j]

      # get the matched sets of SNPs from the background
      snp <- names(d)
      deltas <- abs(snp_freq_sub$dist - d)

      # get matched distance SNPs
      matched_snps <- snp_freq_sub[deltas <= dist_eps, "rsid"]

      return(matched_snps)
    }, mc.cores=threads)

    if(any(unlist(lapply(matched_snps, length)) == 0)) next

    # now do the MAF matching

    # list of MAFs to match, gets reduced each time one MAF is matched
    maf_temp <- mafs
    maf_temp <- maf_temp[!is.na(maf_temp)]
    
    final_snps <- unlist(mclapply(matched_snps, function(d) {
      bg_mafs <- snp_freq_sub[match(d, snp_freq_sub$rsid), "maf"]
      bg_mafs <- bg_mafs[!is.na(bg_mafs)]
      for(j in 1:length(bg_mafs)) {
        deltas <- abs(maf_temp - bg_mafs[j])
        if(any(deltas <= maf_eps)) {
          to_remove <- which.min(deltas)
          maf_temp <<- maf_temp[-to_remove]
          return(d[j])
        }
      }
      return(NA_character_)
    }, mc.cores=threads))

    nas <- sum(is.na(final_snps))

    if (nas == 0) {
      return(list(cpg=icpg, snps=final_snps))
    } else {
      if(nas <= best_match_nas) {
        best_match <- list(cpg=icpg, snps=final_snps)
        best_match_nas <- nas
      }
    }
  }
  warning("Not all SNPs had matches.")
  return(best_match)
}

#' -----------------------------------------------------------------------------
#' Gets a matched background region for the provided CpG and its associated SNPs
#' Special version for use with EPIC result dataframes. Used in the 
#' 'annotation_ocerlap_EPIC.R' script.
#' 
#' We currently sample at most 50,000 background CpGs
#' 
#' -----------------------------------------------------------------------------
get_matched_region_epic <- function(cpg,
                               snps,
                               cosmo_cis,
                               snp_freq_by_chr,
                               cpg_msd,
                               snps_background,
                               cpgs_background,
                               threads = 1) {
  
  require(parallel)
  
  # allow +-1000bp distance mismatch
  dist_eps <- 1000
  
  # allow 0.05 MAF difference
  maf_eps <- 0.05
  
  # distances between CpG and associated SNPs
  a1 <- cosmo_cis[match(snps, cosmo_cis$snp),] %>% pull(snp.pos)
  a2 <- cosmo_cis[match(cpg, cosmo_cis$cpg),] %>% pull(cpg.pos)
  dists <- abs(a1 - a2)
  names(dists) <- snps
  
  # get the mafs of the QTL SNPs
  qtl_snp_chr <- unique(cosmo_cis[match(snps, cosmo_cis$snp),] %>% pull(snp.chr))
  mafs <- snp_freq_by_chr[[qtl_snp_chr]][match(snps,
                                               snp_freq_by_chr[[qtl_snp_chr]]$rsid),] %>% pull(maf)
  names(mafs) <- snps

  matched_cpgs <- rownames(cpg_msd[get.matched.cpgs(cpg, cpg_msd), ])
  matched_cpgs <-
    sample(intersect(cpgs_background$id, matched_cpgs)) # sample cpgs
  
  # get at most 100000 background CpGs
  #max <- 100000
  #if(length(matched_cpgs) > max) {
  #  matched_cpgs <- sample(matched_cpgs, max)
  #}
  
  # get entity position information to match by distance
  cpg_pos <- cpgs_background[matched_cpgs, , drop = F]
 
  # no matched CpGs
  if(nrow(cpg_pos) < 1) {
    return(NA)
  }
  
  best_match <- NA
  best_match_nas <- length(dists)
  
  print(paste0("Checking ", nrow(cpg_pos), " cpgs."))
  
  # check all matched CpGs
  for (i in 1:nrow(cpg_pos)) {
    #print(i)
    icpg <- cpg_pos[i, "id"]
    ichr <- cpg_pos[i, "chr"]
    ipos <- cpg_pos[i, "pos"]
    #  print(ichr)
    # get all SNPs, which are within the max distance
    snp_freq_sub <- subset(snp_freq_by_chr[[ichr]],
                           abs(ipos - pos) <= max(dists) + dist_eps &
                             rsid %in% snps_background$rsid)
    snp_freq_sub$dist <- abs(snp_freq_sub$pos - ipos)
    
    # only proceed if we have as least as many snps as we need within the window
    if(nrow(snp_freq_sub) < length(dists)) next
    
    # for each SNP (named distances) get a matched background SNP such that each
    # distance is reflected once
    matched_snps <- mclapply(1:length(dists), function(j) {
      
      d <- dists[j]
      
      # get the matched sets of SNPs from the background
      snp <- names(d)
      deltas <- abs(snp_freq_sub$dist - d)
      
      # get matched distance SNPs
      matched_snps <- snp_freq_sub[deltas <= dist_eps,] %>% pull(rsid)
      
      return(matched_snps)
    }, mc.cores=threads)
    
    if(any(unlist(lapply(matched_snps, length)) == 0)) next

    # now do the MAF matching
    
    # list of MAFs to match, gets reduced each time one MAF is matched
    maf_temp <- mafs
    maf_temp <- maf_temp[!is.na(maf_temp)]
    
    final_snps <- unlist(mclapply(matched_snps, function(d) {
      bg_mafs <- snp_freq_sub[match(d, snp_freq_sub$rsid),] %>% pull(maf)
      bg_mafs <- bg_mafs[!is.na(bg_mafs)]
      for(j in 1:length(bg_mafs)) {
        deltas <- abs(maf_temp - bg_mafs[j])
        if(any(deltas <= maf_eps)) {
          to_remove <- which.min(deltas)
          maf_temp <<- maf_temp[-to_remove]
          return(d[j])
        }
      }
      return(NA_character_)
    }, mc.cores=threads))
    
    nas <- sum(is.na(final_snps))
    
    if (nas == 0) {
      return(list(cpg=icpg, snps=final_snps))
    } else {
      if(nas <= best_match_nas) {
        best_match <- list(cpg=icpg, snps=final_snps)
        best_match_nas <- nas
      }
    }
    # clear mem cache after a while
    if(i %% 1000 == 0) gc(full=T)
  }
  warning("Not all SNPs had matches.")
  return(best_match)
}
