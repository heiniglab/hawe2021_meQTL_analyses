#' -----------------------------------------------------------------------------
#' Methods for calculating meta-eQTLs based on all LOLIPOP genotype platforms
#' and KORA
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Feb 12 07:17:51 2020
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#' Custom method to generate meta p-values for three studies inside the pmap
#' form purrr. Inverse variance fixed effect meta analysis is employed.
# ------------------------------------------------------------------------------
meta_wrapper <- function(beta1,beta2,beta3,
                         se1,se2,se3) {
  # custom meta analysis (fixed effect, inverse variance)
  #  same results as meta::metagen but LOTS faster (~100 times)
  betas <- c(beta1,beta2,beta3)
  se <- c(se1,se2,se3)
  
  # weights (inv var)
  w <- 1/(se^2)
  sw <- sum(w)
  
  # weighted meta beta
  BETA <- sum(w*betas) / sw
  # weighted meta standard error
  SE <- sqrt(1/sw)
  
  # z-score -> return pvalue
  Z <- BETA/SE
  return(2*pnorm(-abs(Z)))
}

# ------------------------------------------------------------------------------
#' Get meta results from a list of associations containing the KORA EUR and 
#' LOLIPOP EUR/SA results. The dataframe needs to specify three columns with 
#' beta estimates (beta.eur, beta.sa, beta.eurk) and three columns with standard
#' errors (beta_se.eur, beta_se.sa, beta_se.eurk).
#' 
#' @param assocs The dataframe/tibble with the associations
#' @param filter.p The p-value to be used for filtering meta results. 
#' Default taken from the initial enrichment analysis: 5.06e-11
#' 
# ------------------------------------------------------------------------------
get_meta <- function(assocs, filter.p = 5.06e-11) {
 
  if(!require(dplyr) | !require(purrr)) stop("Need dplyr and purrr for this.")
  
  # rowise application of the meta wrapper and instant filter
  meta <- dplyr::mutate(assocs, 
                        meta.p = purrr::pmap_dbl(list(beta.eur, beta.sa, beta.eurk,
                                                      beta_se.eur, beta_se.sa, beta_se.eurk), 
                                                 meta_wrapper)) %>%
    filter(meta.p < filter.p)
  
  return(meta)
}

# ------------------------------------------------------------------------------
#' Adpated version of the functions meta_wrapper and get_meta to return 
#' also the meta beta and meta se values
# ------------------------------------------------------------------------------
meta_wrapper_withBeta <- function(beta1,beta2,beta3,
                                  se1,se2,se3) {
  # custom meta analysis (fixed effect, inverse variance)
  #  same results as meta::metagen but LOTS faster (~100 times)
  betas <- c(beta1,beta2,beta3)
  se <- c(se1,se2,se3)
  
  # weights (inv var)
  w <- 1/(se^2)
  sw <- sum(w)
  
  # weighted meta beta
  BETA <- sum(w*betas) / sw
  # weighted meta standard error
  SE <- sqrt(1/sw)
  
  # z-score -> return pvalue
  Z <- BETA/SE
  tibble::tibble(meta.p=2*pnorm(-abs(Z)),
                 meta.beta=BETA,
                 meta.sd=SE)
}

get_meta_withBeta <- function(assocs, filter.p = 5.06e-11) {
  
  if(!require(dplyr) | !require(purrr)) stop("Need dplyr and purrr for this.")
  
  meta <- purrr::pmap_dfr(list(assocs$beta.eur, assocs$beta.sa, assocs$beta.eurk,
                               assocs$beta_se.eur, assocs$beta_se.sa, assocs$beta_se.eurk), 
                           meta_wrapper_withBeta) %>%
    bind_cols(assocs, .) %>%
    filter(meta.p < filter.p)
  
  return(meta)
}

get_meta_withBeta_woFilter <- function(assocs) {
  
  if(!require(dplyr) | !require(purrr)) stop("Need dplyr and purrr for this.")
  
  meta<-purrr::pmap_dfr(list(assocs$beta.eur, assocs$beta.sa, assocs$beta.eurk,
                            assocs$beta_se.eur, assocs$beta_se.sa, assocs$beta_se.eurk),
                       meta_wrapper_withBeta) %>%
    mutate(cpg=assocs$cpg,
           gene=assocs$gene, .before="meta.p")
  
  return(meta)
}