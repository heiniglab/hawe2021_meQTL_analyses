#' -----------------------------------------------------------------------------
#' Compare the results obtained from the 450k and the EPIC array analyses. uses
#' all outputs from the EPIC results and creates 3 large scatterplots of effects
#' for 450k vs EPIC
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Jan  09 14:32:38 2019
#' -----------------------------------------------------------------------------

library(tidyverse)
source("scripts/lib.R")
library(ggplot2)
library(cowplot)
library(parallel)
theme_set(theme_cowplot())

# load data
disco_450k <- load_st1_pairs(fst1="data/benni.lehne/st1_v2.txt")
epic_in_st1 <- read_tsv("results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/all_associations_in_st1.tsv")
common_probes <- read_tsv("results/current/common_probes.txt", col_names=F) %>% pull(X2)

# filter for common probes only
disco_450k_filtered <- disco_450k %>%
  filter(cpg %in% common_probes)
epic_filtered <- epic_in_st1 %>%
  filter(cpg %in% disco_450k_filtered$cpg)

print("Joining results.")
# for each associations on the LHS, we add corresponding associations from the RHS
# if no match is possible, fields will be populated with NA for the RHS
joined <- left_join(disco_450k_filtered, epic_filtered,
                    by = c("snp" = "SNP", "cpg" = "cpg"),
                    suffix = c(".450k",".epic")) %>%
  mutate(is_epic = !is.na(beta) & `p-value` < 0.05)
  mutate(has_epic = !is.na(beta))

# create the scatter plot
gp <- ggplot(joined, aes(y=beta.eur.discovery, x=beta)) + 
  geom_bin2d(bins=1000) + 
  scale_fill_gradient(low="darkblue", high="orange") +
  facet_wrap(.~ category, ncol=2) +
  background_grid(major = "xy") +
  labs(x="EPIC beta",y="DISCOVERY beta", title = "Scatterplot (binned) of effect sizes.")

# ------------------------------------------------------------------------------
print("Saving plot.")
# ------------------------------------------------------------------------------
save_plot(filename = "results/current/scatterplot.pdf", gp, ncol=2, nrow=2)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
