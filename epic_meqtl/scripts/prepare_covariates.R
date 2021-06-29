#' -----------------------------------------------------------------------------
#' Prepare covariates for the association analysis. 
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Dec  6 14:32:38 2019
#' -----------------------------------------------------------------------------


log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# inputs
fidmap <- snakemake@input$idmap
fgeno_indiv <- snakemake@input$geno_indiv
fhouseman <- snakemake@input$houseman
finput_covariates <- snakemake@input$covariates
ftechnical_variables <- snakemake@input$technical_variables
fplatebatch <- snakemake@input$platebatch
findiv_ignore <- snakemake@input$indiv_ignore

# outputs
fgeno_id_ordering <- snakemake@output$geno_id_ordering
fgeno_ids <- snakemake@output$geno_ids
fmeth_ids <- snakemake@output$meth_ids
foutput_covariates <- snakemake@output$covariates

# ------------------------------------------------------------------------------
print("Load ID map and covariates, using only individuals with meth/geno data.")
# ------------------------------------------------------------------------------
# methylation individuals to be ignored
indiv_ignore <- read_tsv(findiv_ignore, col_names=F, 
                         col_types=cols(col_character())) %>% pull(X1)

# id mapping
IDMAP <- read_delim(fidmap,
                    delim=";",
                    col_types = cols(.default=col_character())) %>%
    rename(meth_id = u3n_meth850K_ff4, geno_id = lg_dnaAxiom_s4f4) %>%
    select(meth_id, geno_id) %>%
    drop_na() %>%
    filter(!meth_id %in% indiv_ignore)

# genotype individuals
geno_indiv <- read_tsv(fgeno_indiv, col_names=F,
                       col_types=cols(.default=col_character())) %>%
    rename(geno_id = X1)

# load houseman data and remove Gran
houseman <- read_delim(fhouseman,
                       col_types=cols(zz_nr_ff4_methyl = col_character()),
                       delim = ";") %>%
    rename(meth_id=zz_nr_ff4_methyl) %>%
    select(-Gran)

# covariates
covars <- read_delim(finput_covariates,
                     delim=",",
                     col_types=cols(zz_nr_ff4_methyl=col_character(),
                                    u3l_mdif=col_character())) %>%
    rename(meth_id=zz_nr_ff4_methyl, age = u3talteru,
           sex = u3csex, bmi = u3tbmi, wbc = u3l_wbc) %>%
    select(meth_id,age,sex,bmi,wbc)

# batch variable
load(fplatebatch)
batch <- PlateChip %>% as_tibble() %>%
  rename(meth_id = ZZ_nr) %>%
  mutate(meth_id = as.character(meth_id)) %>%
  select(meth_id, Plate, Batch)

# control probe PCs
load(ftechnical_variables)
meth_pcs <- data.frame(pcs)[,1:20]
meth_pcs$meth_id <- rownames(meth_pcs)

# combine IDMAP and covariates
IDMAP_covars <- IDMAP %>%
  inner_join(y=geno_indiv, by=c("geno_id" = "geno_id")) %>%
  inner_join(y=houseman, by=c("meth_id" = "meth_id")) %>%
  inner_join(y=covars, by=c("meth_id" = "meth_id")) %>%
#  inner_join(y=meth_pcs, by=c("meth_id" = "meth_id")) %>%
  inner_join(y=batch, by=c("meth_id" = "meth_id"))

print(paste0("Number of individuals with covariates: ", nrow(IDMAP_covars)))

# ------------------------------------------------------------------------------
print("All done, saving.")
# ------------------------------------------------------------------------------

# report (ordered) list of geno ids
write_tsv(IDMAP_covars %>% select(geno_id), fgeno_ids,
          col_names = FALSE)

# report actual ordering for genotype file processing
ordering <- match(IDMAP_covars %>% pull(geno_id), geno_indiv %>% pull(geno_id))
write_tsv(enframe(name=NULL, ordering), fgeno_id_ordering,
          col_names = FALSE)

# report all covariates
write_tsv(IDMAP_covars %>% select(meth_id), fmeth_ids,
          col_names = FALSE)

# report all covariates
out <- IDMAP_covars %>% select(-meth_id, -geno_id)

# create model mastrix for matrixEQTL, remove intercept definition
out <- model.matrix(~age+sex+bmi+wbc+CD8T+CD4T+NK+Bcell+Mono+Plate, out)
out <- out[, !grepl("Intercept", colnames(out))]
out <- t(as.data.frame(out))

colnames(out) <- IDMAP_covars$geno_id

write.table(as.data.frame(out), foutput_covariates,
            row.names=T, col.names=NA, quote=F, sep="\t")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
