# ------------------------------------------------------------------------------
# Extract from the full set of EPIC associations only those which match
# directly (SNP id) to one of the sentinel hotspots.
# ------------------------------------------------------------------------------


library(tidyverse)
source("scripts/lib.R")

epic_assocs <- readRDS("results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/all_associations_discovery_annotated.rds")
sentinels <- read_tsv("data/sentinels.txt", col_names=F) %>% pull(X1)

filtered <- filter(epic_assocs, SNP %in% sentinels & category == "trans")

print("Number of new associations:")
print(nrow(filtered))

# filter for CpGs with TFBS 
tfbs_annot <- load_tfbs_annot(unique(filtered$cpg), 100)
cpg_has_tfbs <- rowSums(tfbs_annot)>1

filtered <- filtered %>% mutate(cpg_has_tfbs = cpg_has_tfbs[cpg])

# add with information whether cpgs are epic specific
array_annot <- read_tsv("results/current/cpg_annotation_both_arrays.tsv")

filtered_annot <- left_join(filtered, array_annot, by=c("cpg"="id")) %>% 
  dplyr::rename(on_epic = is_epic, on_450k = is_450k, on_both = is_both, 
  cpg_mean_450k = mean.450k, cpg_sd_450k = sd.450k,
  cpg_mean_epic = mean.epic, cpg_sd_epic = sd.epic)

write_tsv(filtered_annot, 
  "results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/sentinel_associations_discovery_annotated.tsv")


