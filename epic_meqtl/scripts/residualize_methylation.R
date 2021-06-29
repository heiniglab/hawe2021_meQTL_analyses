#' -----------------------------------------------------------------------------
#' get residualized data matrix for EPIC methylation data.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Jan  28 14:32:38 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(parallel)
library(GenomicRanges)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
fmeth <- snakemake@input$methylation
fhouseman <- snakemake@input$houseman
finput_covariates <- snakemake@input$covariates
fidmap <- snakemake@input$idmap
fplatebatch <- snakemake@input$platebatch
findiv_ignore <- snakemake@input$indiv_ignore

# outputs
fresult <- snakemake@output$result
fresult_no_houseman <- snakemake@output$result_no_houseman

# ------------------------------------------------------------------------------
print("Loading data..")
# ------------------------------------------------------------------------------
# methylation data
load(fmeth)

# methylation individuals to be ignored
indiv_ignore <- read_tsv(findiv_ignore, col_names=F,
                         col_types=cols(col_character())) %>% pull(X1)

# id mapping
IDMAP <- read_delim(fidmap,
                    delim=";",
                    col_types = cols(.default=col_character())) %>%
    dplyr::rename(meth_id = u3n_meth850K_ff4, geno_id = lg_dnaAxiom_s4f4) %>%
    dplyr::select(meth_id) %>%
    drop_na() %>%
    filter(!meth_id %in% indiv_ignore)

# load houseman data and remove Gran
houseman <- read_delim(fhouseman,
                       col_types=cols(zz_nr_ff4_methyl = col_character()),
                       delim = ";") %>%
    dplyr::rename(meth_id=zz_nr_ff4_methyl) %>%
    dplyr::select(-Gran)

# covariates
covars <- read_delim(finput_covariates,
                     delim=",",
                     col_types=cols(zz_nr_ff4_methyl=col_character(),
                                    u3l_mdif=col_character(),
                                    u3l_plaz=col_character(),
				    u3l_igr=col_character())) %>%
    dplyr::rename(meth_id=zz_nr_ff4_methyl, age = u3talteru,
           sex = u3csex, bmi = u3tbmi, wbc = u3l_wbc) %>%
    dplyr::select(meth_id,age,sex,bmi,wbc)

# batch variable
load(fplatebatch)
batch <- PlateChip %>% as_tibble() %>%
  dplyr::rename(meth_id = ZZ_nr) %>%
  mutate(meth_id = as.character(meth_id)) %>%
  dplyr::select(meth_id, Plate, Batch)

# combine IDMAP and covariates
IDMAP_covars <- IDMAP %>%
  inner_join(y=houseman, by=c("meth_id" = "meth_id")) %>%
  inner_join(y=covars, by=c("meth_id" = "meth_id")) %>%
  inner_join(y=batch, by=c("meth_id" = "meth_id"))

print(paste0("Number of individuals with covariates: ", nrow(IDMAP_covars)))

# report all covariates
mmatrix <- IDMAP_covars %>% select(-meth_id)
# create model matrix for matrixEQTL, remove intercept definition
mmatrix <- model.matrix(~age+sex+bmi+wbc+CD8T+CD4T+NK+Bcell+Mono+Plate, mmatrix)
rownames(mmatrix) <- IDMAP_covars$meth_id

# and the one without white blood-cell counts and houseman estimates:
mmatrix_no_hm <- IDMAP_covars %>% select(-meth_id)
mmatrix_no_hm <- model.matrix(~age+sex+bmi+Plate, mmatrix_no_hm)
rownames(mmatrix_no_hm) <- IDMAP_covars$meth_id

# subset methylation data to our individuals
beta <- t(beta[,IDMAP_covars %>% pull(meth_id)])

# ------------------------------------------------------------------------------
print("Getting residuals.")
# ------------------------------------------------------------------------------
beta_resid <- apply(beta, 2, function(cpg) {
  resid(lm(cpg ~ mmatrix, na.action=na.exclude))
})
print("Saving.")
saveRDS(beta_resid, fresult)

rm(beta_resid)
gc()


print("no houseman estimates and wbc...")
beta_resid <- apply(beta, 2, function(cpg) {
  resid(lm(cpg ~ mmatrix_no_hm, na.action=na.exclude))
})
print("Saving.")
saveRDS(beta_resid, fresult_no_houseman)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
