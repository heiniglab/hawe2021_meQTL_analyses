#' -----------------------------------------------------------------------------
#' Library for the long-range/trans Hi-C and TAD enrichment analyses
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Dec 16 11:48:54 2019
#' -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Load the provided TAD regions definitions into a GRanges object
#'
#' @param ftads The file from which to load the TADs
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
load_tad_regions <- function(ftads) {
  require(GenomicRanges)
  require(readr)
  
  tads <- read_tsv(ftads, col_names=T)
  rtads <- with(tads, GRanges(paste0("chr", chr),
                              IRanges(start,end)))
  
  return(rtads)
}

# ------------------------------------------------------------------------------
#' Load the HiC regions into two GRanges objects (i.e. target/source, etc.)
#'
#' Will return a list of two GRanges objects with names 'r1' and 'r2'
#' 
#' @param fhic Path to file containing the HiC contacts.
#' @param fpchic Path to file containing the PCHiC contacts.
#' @param extension Amount of BP by which contact regions shouls be extendeded in
#' both directions
#' @param resolution Resolution of the HiC data, not relevant for PCHIC data
#' @param is.pchic Flag whether we get promoter capture HiC data which has to
#' be processes a little differentily
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
load_hic_regions <- function(fhic = NULL, fpchic = NULL,
                             extension = 0,
                             resolution = 0) {

  require(GenomicRanges)
  require(readr)
  
  if(is.null(fhic) & is.null(fpchic)) stop("No data given")
  if(!is.null(fhic) & !is.null(fpchic)) stop("Only one data file can be specified")
  
  if(!is.null(fhic)) {
    hic <- read_tsv(fhic, col_names=F, col_types = 
                      cols(.default = col_double(), 
                           X2 = col_character(), X4 = col_character()))
  
    colnames(hic) <- c("count", "chr1", "start1", "chr2", "start2")
    r1 <- with(hic, GRanges(paste0("chr", chr1),
                              IRanges(start1 - extension,
                                      start1 + resolution + extension)))
    r2 <- with(hic, GRanges(paste0("chr", chr2),
                               IRanges(start2 - extension,
                                       start2 + resolution + extension)))
  } else {
    # pchic data
    pchic <- read_tsv(fpchic, col_types = cols(oeChr = col_character(), baitChr = col_character())) %>%
      filter(baitChr == oeChr)
    
    r1 <- with(pchic, GRanges(paste0("chr", baitChr), IRanges(baitStart-extension, 
                                                              baitEnd+extension)))
    r2 <- with(pchic, GRanges(paste0("chr", baitChr), IRanges(oeStart-extension, 
                                                              oeEnd+extension)))
  }
  
  return(list(r1=r1, r2=r2))
}


# ------------------------------------------------------------------------------
#' Method to quickly check whether two entities are in HiC contact regions.
#' e1 and e2 as well as source and target need to be 'parallel' ranges, i.e.
#' for each i, e1[i] and e2[i] correspond to the same QTL pair, and for each j
#' source[j] and target[j] correspond to the same contact-
#'
#' @param e1 The first entity to be checked (usually a SNP)
#' @param e2 The second entity to be checked (usually a CpG)
#' @param source The contact HiC regions' source regions (e.g. bait)
#' @param target The contact HiC regions' target regions (e.g. 'other')
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
any_in_same_contact <- function(e1, e2, source, target) {
  require(dplyr)
  
  # get all overlaps
  ov1 <- as(findOverlaps(source, e1), "data.frame")
  ov2 <- as(findOverlaps(target, e2), "data.frame")
  
  # inner join by default checks all columns with the same name for matches
  matches <- inner_join(ov1, ov2, by = c("queryHits" = "queryHits", 
                                         "subjectHits" = "subjectHits"))
  
  # extract query and subject hits for quick check
  i <- NA
  if(nrow(matches) > 0) {
    # return first hit
    hit <- TRUE
    i <- matches[1, "subjectHits"]
  } else {
    # check other way around
    ov1 <- as(findOverlaps(source, e2), "data.frame")
    ov2 <- as(findOverlaps(target, e1), "data.frame")
    
    # inner join by default checks all columns with the same name for matches
    matches <- inner_join(ov1, ov2, by = c("queryHits" = "queryHits", 
                                           "subjectHits" = "subjectHits"))
    
    if(nrow(matches) > 0) {
      # return first hit
      hit <- TRUE
      i <- matches[1, "subjectHits"]
    } else {# return random pair
      hit <- FALSE
      i <- 1
    }
  }
  names(hit) <- paste0(names(e1)[i], "_", names(e2)[i])
  return(hit)
}

# ------------------------------------------------------------------------------
#' Method to quickly check whether pairs of entities are within the same TAD
#'
#' e1 and e2 need to be parallel vectors. 
#' Compare also methods 'any_in_same_contact'
#' 
#' @param e1 The first entity list to be checked (e.g. list of SNPs)
#' @param e2 The second entity list to be checked (e.g. list of CpGs)
#' @param tads The TAD regions as GRanges
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
any_in_same_tad <- function(e1, e2, tads) {
  # get all overlaps
  ov1 <- as(findOverlaps(tads, e1), "data.frame")
  ov2 <- as(findOverlaps(tads, e2), "data.frame")
  
  # inner join by default checks all columns with the same name for matches
  matches <- inner_join(ov1, ov2, by = c("queryHits" = "queryHits", 
                                         "subjectHits" = "subjectHits"))
  
  # extract query and subject hits for quick check
  i <- NA
  if(nrow(matches) > 0) {
    # return first hit
    hit <- TRUE
    i <- matches[1, "subjectHits"]
  } else {
    hit <- FALSE
    i <- 1
  }
  names(hit) <- paste0(names(e1)[i], "_", names(e2)[i])
  return(hit)
}


# ------------------------------------------------------------------------------
#' Extract neat enrichment results from the generated enrichment list
#' 
#' @param result The list of enrichment results
#' @param index Name of 'foreground' enrichment elements for the items in the 
#' enrichment list
#' @param index_bg Name of 'background' enrichment elements for the items in the 
#' enrichment list
#' 
# ------------------------------------------------------------------------------
extract_enrichment <- function(result, index, index_bg) {
  
  # get hits for index
  index_hits <- sapply(result, "[[", index)
  index_summary <- table(index_hits)

  # get the background summary
  bg_hits <- unlist(lapply(result, "[[", index_bg))
  bg_summary <- table(bg_hits, useNA="no")
  
  # report number of NAs
  print(paste0("Missing backgrounds:", length(which(is.na(bg_hits)))))
  
  # build contingency table
  cont <- get_cont_table(index_hits, bg_hits)
  
  # calculate fisher test and get p-value
  pv <- fisher.test(cont,
                    alternative="g")$p.value
  
  print("Results:")
  print(cont)
  print(paste0("p-value:", pv))
  
  # create final result list
  result <- list(cont,
                 pv,
                 index_hits,
                 index_summary,
                 bg_hits,
                 bg_summary)
  names(result)<- c("contingency_table", 
                    "pvalue",
                    "index_hits", "index_summary",
                    "bg_hits", "bg_summary")
  result
}

# ------------------------------------------------------------------------------
# Get matched SNPs (MAF) for a specific rsid
# ------------------------------------------------------------------------------
get_matched_snps <- function(snpid, snp_annotation, snp_ranges, 
                             bg_snps) {
  
  matched_snps <- get.matched.snps(snpid, snp_annotation, af_column = "maf")
  matched_snps <- intersect(matched_snps, names(bg_snps))
  matched_snps <- sample(matched_snps,
                         min(50000, length(matched_snps)))
  
  matched_snps_ranges <- bg_snps[matched_snps]

  return(matched_snps_ranges)
  
}

# ------------------------------------------------------------------------------
# Enrich longrange eQTLs in TAD and HiC regions.
#'
#' @param pairing List of lists. Contains for each 'sentinel', a list of SNPs
#' and Cpgs which can be mapped to the sentinel. It is possible to work without
#' the definition of sentinels, by simply supplying 1:1 lists (1 SNP : 1 CpG)
#' @param snp_ranges For all candidate SNPs to be tested, the corresponding 
#' genomic information as a GRanges object. Needs to have names matching the 
#' candidate SNPs.
#' @param cpg_ranges As 'snp_ranges', but for the CpGs.
#' @param tad_regions The TAD regions in which to enrich the pairs
#' @param hic_source The HIC contact 'source' regions (first region of contact)
#' @param hic_target The HIC contact 'target' regions (second region of contact)
#' @param snp_annotation data.frame of SNP info (containing MAF and positions of
#' all BG and QTL SNPs)
#' @param snp_annotation data.frame of CpG info (containing mean betas and 
#' position info of all BG and QTL CpGs)
#' @param bg_cpgs GRanges object containing all background CpGs. Names must be
#' CpG ids
#' @param bg_snps GRanges object containing all background SNPs. Names must be
#' SNP rsids
#' @param dist_tolerance Tolerance for distance matching of background SNP/CpG
#' pairs
#' @param threads Number of threads to be used during enrichment analysis, 
#' defaults to 1.
#' 
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
enrich_longrange <- function(pairing,
                         snp_ranges,
                         cpg_ranges,
                         tad_regions,
                         hic_source,
                         hic_target,
                         snp_annotation,
                         cpg_annotation,
                         bg_cpgs,
                         bg_snps,
                         dist_tolerance,
                         threads = 1) {
  
  require(tidyr) # crossing()
  
  bg_cpg_ids <- names(bg_cpgs)
  
  # iterate over all pairing
  temp <- mclapply(1:length(pairing), function(i) {
    
  #temp <- mclapply(21:100, function(i) {
    # get the ids for this instance
    snp_ids <- pairing[[i]]$snps
    cpg_ids <- pairing[[i]]$cpgs
    
    # we want to check all combinations of SNPs and CpGs
    # ( e.g. for a single sentinel pair mapped back to original pairs)
    cross <- crossing(snp_ids, cpg_ids)
    
    # get entities
    snps <- snp_ranges[cross$snp_ids]
    cpgs <- cpg_ranges[cross$cpg_ids]
    
    # check overlap of meQTL with TAD and HiC
    tad_overlap <- suppressWarnings(
      any_in_same_tad(snps, cpgs, tad_regions)
    )
    hic_overlap <- suppressWarnings(
      any_in_same_contact(snps, cpgs, hic_source, hic_target)
    )
    
    # the overlap results
    tad_overlap_bg <- NULL
    hic_overlap_bg <- NULL
    
    # matched entities
    msnp_contact <- NULL
    mcpg_contact <- NULL
    msnp_tad <- NULL
    mcpg_tad <- NULL
    
    # Get background (non-meQTL) pairs individually for all combinations
    for(j in 1:nrow(cross)) {
      
      snp <- snps[j]
      cpg <- cpgs[j]
      
      # get all matched SNPs
      matched_snps_ranges <- get_matched_snps(names(snp), snp_annotation,
                                              snp_ranges, bg_snps)
      
      # get distance to original snp and sample all cpgs available with same
      # distance to the sampled SNP
      d <- distance(snp,cpg)
      # define lower and upper boundries for distance matching
      dl <- d - dist_tolerance
      du <- d + dist_tolerance
    
      # get CpGs matched by beta, apply distance matching below
      matched_cpgs <- get.matched.cpgs(names(cpg), cpg_annotation)
      matched_cpgs <- rownames(cpg_annotation)[matched_cpgs]
      matched_cpgs <- intersect(matched_cpgs, bg_cpg_ids)
      
      # iterate over all SNPs and return the matching pairs
      for(k in 1:length(matched_snps_ranges)) {
        
        msnp <- matched_snps_ranges[k]
        
        # get distance matched CpGs
        scpgs <- get_distance_matched_cpgs(msnp, bg_cpgs,
                                           dl, du, matched_cpgs)
        scpgs_promoter <- scpgs[scpgs$in_promoter]
        
        # we got some CpGs!
        if(length(scpgs)>0) {
          
          if(is.null(tad_overlap_bg)) {
            # in any case, get the TAD overlap for one of the CpGs
            mcpg_tad <- sample(scpgs, 1)
            msnp_tad <- msnp
            tad_overlap_bg <- any_in_same_tad(msnp_tad, mcpg_tad, tad_regions)
          }
          
          # for HiC enrichment, also check that at least one entity is in a 
          # promoter region
          
          # snp in promoter, just sample any CpG
          if((msnp$in_promoter) & runif(1)<0.5) {
            msnp_contact <- msnp
            mcpg_contact <- sample(scpgs, 1)
            hic_overlap_bg <- any_in_same_contact(msnp_contact,
                                                  mcpg_contact,
                                                  hic_source, hic_target)
            break
          } else {
            # matched promoter, stop now
            if((length(scpgs_promoter) > 0)) {
              msnp_contact <- msnp
              mcpg_contact <- sample(scpgs_promoter, 1)
              hic_overlap_bg <- any_in_same_contact(msnp_contact,
                                                    mcpg_contact,
                                                    hic_source, hic_target)
              break
            }
          }
        }
      }
      
      # we exit if we matched both categories
      if(!is.null(tad_overlap_bg) & !is.null(hic_overlap_bg)) {
        break
      }
    }
    
    # --------------------------------------------------------------------------
    # remember names and finalize
    # --------------------------------------------------------------------------
    if(!is.null(tad_overlap_bg)) {
      names(tad_overlap_bg) <- paste0(names(msnp_tad), "_",
                                      names(mcpg_tad), "_",
                                      distance(msnp_tad, mcpg_tad))
    } else {
      tad_overlap_bg <- NA
    }
    if(!is.null(hic_overlap_bg)) {
      names(hic_overlap_bg) <- paste0(names(msnp_contact), "_",
                                      names(mcpg_contact), "_",
                                      distance(msnp_contact, mcpg_contact))
    } else {
      hic_overlap_bg <- NA
    }
    
    # create result list
    list(hit_in_tad = tad_overlap,
         hit_in_contact = hic_overlap,
         bg_in_tad = tad_overlap_bg,
         bg_in_contact = hic_overlap_bg)
    
  }, mc.cores=threads)
  
  # ----------------------------------------------------------------------------
  print("Extract enrichments for QTLs and background.")
  # ----------------------------------------------------------------------------
  enrichment_tad <- extract_enrichment(temp, 
                                       "hit_in_tad", "bg_in_tad") 
  enrichment_hic <- extract_enrichment(temp, 
                                       "hit_in_contact", "bg_in_contact") 
  
  # build and return all results
  list(tad_enrich = enrichment_tad, 
       hic_enrich = enrichment_hic)
}

# ------------------------------------------------------------------------------
#' Method DESC
#'
#' @param pairing List of lists. Contains for each 'sentinel', a list of SNPs
#' and Cpgs which can be mapped to the sentinel. It is possible to work without
#' the definition of sentinels, by simply supplying 1:1 lists (1 SNP : 1 CpG)
#' @param snp_ranges For all candidate SNPs to be tested, the corresponding 
#' genomic information as a GRanges object. Needs to have names matching the 
#' candidate SNPs.
#' @param cpg_ranges As 'snp_ranges', but for the CpGs.
#' @param hic_source The HIC contact 'source' regions (first region of contact)
#' @param hic_target The HIC contact 'target' regions (second region of contact)
#' @param snp_annotation data.frame of SNP info (containing MAF and positions of
#' all BG and QTL SNPs)
#' @param snp_annotation data.frame of CpG info (containing mean betas and 
#' position info of all BG and QTL CpGs)
#' @param bg_cpgs GRanges object containing all background CpGs. Names must be
#' CpG ids
#' @param bg_snps GRanges object containing all background SNPs. Names must be
#' SNP rsids
#' @param threads Number of threads to be used during enrichment analysis, 
#' defaults to 1.
#' 
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
enrich_trans <- function(pairing,
                         snp_ranges,
                         cpg_ranges,
                         hic_source,
                         hic_target,
                         snp_annotation,
                         cpg_annotation,
                         bg_cpgs,
                         bg_snps,
                         threads = 1) {
  
  require(tidyr) # crossing()
  
  enrichment <- mclapply(1:length(pairing), function(i) {
  #enrichment <- mclapply(600:800, function(i) {
    # get the ids for this instance
    snp_ids <- pairing[[i]]$snps
    cpg_ids <- pairing[[i]]$cpgs
    
    # we want to check all combinations of SNPs and CpGs
    # ( e.g. for a single sentinel pair mapped back to original pairs)
    cross <- crossing(snp_ids, cpg_ids)

    # get entities
    snps <- snp_ranges[cross$snp_ids]
    cpgs <- cpg_ranges[cross$cpg_ids]
    
    # check overlap of meQTL with HiC
    meqtl_overlap <- suppressWarnings(
      any_in_same_contact(snps, cpgs, hic_source, hic_target)
      )
    
    # Get background (non-meQTL) pairs individually for all combinations
    for(j in 1:nrow(cross)) {
    
      snp <- snps[j]
      cpg <- cpgs[j]
      
      # get CpGs matched by beta, apply distance matching below
      matched_cpgs <- get.matched.cpgs(names(cpg), cpg_annotation)
      matched_cpgs <- rownames(cpg_annotation)[matched_cpgs]
      matched_cpgs <- intersect(matched_cpgs, names(bg_cpgs))
      
      # sample background SNPs matched by AF
      matched_snps <- get.matched.snps(names(snp), snp_annotation, 
                                       af_column = "maf")
      matched_snps <- intersect(matched_snps, names(bg_snps))
      matched_snps <- sample(matched_snps,
                             min(50000, length(matched_snps)))
      matched_snps_ranges <- bg_snps[matched_snps]
      
      # the overlap results
      bg_overlap <- NA
      
      # matched contact entities
      msnp <- NULL
      mcpg <- NULL
      
      # no initial matches -> nothing to do for this pair
      if(length(matched_snps) < 1) next
      if(length(matched_cpgs) < 1) next
      
      # iterate over all SNPs and return the matching pairs
      for(k in 1:length(matched_snps_ranges)) {
        
        msnp <- matched_snps_ranges[k]
        
        # get 'distance matched' CpGs (ensures different chromosome to snp)
        scpgs <- get_distance_matched_cpgs(msnp, bg_cpgs,
                                           -1, -1, matched_cpgs, 
                                           is.trans=T)
        # we got some CpGs!
        if(length(scpgs)>0) {
          # sample one of the matched CpGs, SNPs got already sampled at the 
          # beginning
          mcpg <- sample(scpgs, 1)
          bg_overlap <- suppressWarnings(
            any_in_same_contact(msnp,
                                mcpg,
                                hic_source, 
                                hic_target)
            )
          break
        }
      }
      
      # set names and check whether to set NA result
      if(!is.na(bg_overlap)) {
        names(bg_overlap) <- paste0(names(msnp), "_",
                                    names(mcpg), "_")
      }
    }
    
    list(meqtl_overlap = meqtl_overlap,
         bg_overlap = bg_overlap)
    
  }, mc.cores=threads)
    
  # ----------------------------------------------------------------------------
  print("Extracting enrichment, finalizing.")
  # ----------------------------------------------------------------------------
  result <- extract_enrichment(enrichment, "meqtl_overlap", "bg_overlap")
  
  return(result)
}
