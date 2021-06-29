# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(reshape2)
library(parallel)
theme_set(theme_cowplot())
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
meth <- readRDS("results/current/methylation_residualized.rds")
sentinels <- read_tsv("data/benni.lehne/all_r02_110417_combined.tsv")
cpg_annot <- read_tsv("results/current/cpg_annotation_both_arrays.tsv")

# load cpg ranges for distance matching
cpg_ranges <- get_850k_cpg_ranges()

# get unique sentinel cpg ranges
sentinel_cpgs <- sentinels %>% pull(cpg.sentinel) %>% unique
sentinel_cpg_ranges <- cpg_ranges[sentinel_cpgs[sentinel_cpgs %in% names(cpg_ranges) &
                                                sentinel_cpgs %in% colnames(meth)]]

# get epic cpg ranges
epic_cpgs <- filter(cpg_annot, is_epic & !is_450k) %>% pull(id)
epic_cpg_ranges <- cpg_ranges[epic_cpgs]

# ------------------------------------------------------------------------------
# for each EPIC cpg, get the sentinel CPG nearest to it
# ------------------------------------------------------------------------------
nearest_to_epic <- sentinel_cpg_ranges[nearest(epic_cpg_ranges, 
                                               sentinel_cpg_ranges)]

# ------------------------------------------------------------------------------
# get all r^2s for the nearest entities
# ------------------------------------------------------------------------------
check_cpg_r2 <- function(c1, c2) {
  meth_r2 <- cor(meth[,c1], meth[,c2],
				 use="complete")^2
  return(meth_r2)
}

threads <- 12

r2s <- unlist(mclapply(1:length(epic_cpg_ranges), function(i) { 
  check_cpg_r2(names(epic_cpg_ranges)[i], names(nearest_to_epic)[i])
}, mc.cores = threads))

table(r2s>0.2)
#
# FALSE   TRUE
#374467   7138

# ------------------------------------------------------------------------------
# Check distances
# ------------------------------------------------------------------------------
dist2nearest <- distanceToNearest(epic_cpg_ranges, sentinel_cpg_ranges)
table(mcols(dist2nearest)$distance < 1e6)
#
# FALSE   TRUE
#   586 381019

# ------------------------------------------------------------------------------
# create a smoothscatter plot for r2s VS distance
# probably not necessary to reorder, but do anyways
# ------------------------------------------------------------------------------
epic_cpg_ranges <- epic_cpg_ranges[queryHits(dist2nearest)]
r2s <- r2s[queryHits(dist2nearest)]
dists <- mcols(dist2nearest)[queryHits(dist2nearest),"distance"]

toplot <- tibble(r2=r2s, distance=dists)

write_tsv(toplot, "results/current/EPIC_to_sentinel_CpG_r2s.tsv")

# final plot used in response
gp2 <- toplot %>%
  ggplot(aes(x=log10(distance+1), y=r2), toplot) +
  scale_fill_gradient(trans="log", breaks=c(1,7,55,400,3000)) + 
  geom_bin2d(bins=200)

# Full r2s were calculated in separate script!
r2s_no_hm <- read_tsv("results/current/EPIC_r2s_with_distances_no_housemand_and_wbc.tsv")
r2s <- read_tsv("results/current/EPIC_r2s_with_distances.tsv")
zoomed <- r2s %>% ggplot(aes(x=distance, y=r2)) + geom_bin2d(bins=10000) + 
  coord_cartesian(xlim=c(1,6e3)) + scale_fill_gradient(trans="log", breaks=c(1,7,55,400,3500))
zoomed_no_hm <- r2s_no_hm %>% ggplot(aes(x=distance, y=r2)) + geom_bin2d(bins=10000) + 
  coord_cartesian(xlim=c(1,6e3)) + scale_fill_gradient(trans="log", breaks=c(1,7,55,400,3500))

final <- plot_grid(gp2+labs(title="EPIC vs closest sentinel"), 
                   zoomed+labs(title="EPIC vs EPIC, r-squared"),
                   zoomed_no_hm+labs(title="EPIC vs EPIC, r-squared (no Houseman)"), ncol=2, labels=c("A", "B1", "B2"))
save_plot("results/current/epic_summary.pdf", final, ncol=2, nrow=2)


# cpgs with r2>0.2 and distance < 50 -> check for binding sites?
dist_cutoff <- 50

cpgs_in_window <- names(epic_cpg_ranges[dists<=dist_cutoff & r2s > 0.2])
cpgs_outside_window <- names(epic_cpg_ranges[dists>dist_cutoff & r2s > 0.2])

print("Number of CpGs within 100bp window:")
print(length(cpgs_in_window))

print("Number of CpGs ouside 100bp window:")
print(length(cpgs_outside_window))

cpgs_to_check <- cpgs_outside_window

write_tsv(tibble(cpgs=cpgs_to_check), paste0("results/current/cpgs_to_check_with_tfbs_", 
                                             dist_cutoff, ".txt"), 
          col_names=F)

# get tfbs context
tfbs_annot <- load_tfbs_annot(cpgs_to_check, 100)

ncpgs_with_tfbs <- sum(apply(tfbs_annot, 1, any))
print("Number of CpGs outside window with TFBS:")
print(ncpgs_with_tfbs)

cpgs_with_tfbs <- rownames(tfbs_annot[apply(tfbs_annot, 1, any),])
write_tsv(tibble(cpgs=cpgs_with_tfbs), col_names=F, path="cpgs_with_tfbs.txt")

# get trans-meQTL from the above list of CpGs with r2>0.2 (with nearest sentinel) CpGs
# bash:
# awk '{if($5+0 < 1e-14) print;}' results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/trans_associat$
#  fgrep -w -f results/current/cpgs_to_check_with_tfbs.txt | cut -f 2 | sort | uniq > trans_cpgs_to_check.txt
# or (for CpGs with TFBS):
#  fgrep -w -f cpgs_with_tfbs.txt | cut -f 2 | sort | uniq > trans_cpgs_with_tfbs.txt

# further process the interesting trans meQTL CpGs with tfbs
trans_cpgs <- read_tsv("trans_cpgs_with_tfbs.txt", col_names=F) %>% pull(X1)
map <- tibble(epic=names(epic_cpg_ranges), sentinel=names(sentinel_cpg_ranges[subjectHits(dist2nearest)])) %>%
   filter(epic %in% trans_cpgs)

write_tsv(path="trans_cpgs_with_tfbs_mapping.tsv", map)

# check the missing CpGs for their categories
missing <- setdiff(sentinel_cpgs, names(sentinel_cpg_ranges))
sentinels %>% filter(cpg.sentinel %in% missing) %>% select(cpg.sentinel, category) %>% distinct %>% pull(category) %>% table


# check all nearest distances and the respective correlations
meth <- readRDS("results/current/methylation_residualized.rds")
cpg_ranges_sub <- cpg_ranges[names(cpg_ranges) %in% colnames(meth) & grepl("^cg", names(cpg_ranges)]
nearest_ranges <- cpg_ranges_sub[nearest(cpg_ranges_sub)]

# create data frame
toplot <- tibble(cpg1 = names(cpg_ranges_sub), cpg2 = names(nearest_ranges), distance = distance(cpg_ranges_sub, nearest_ranges))

# get r2s for each of the pairs, this takes a few hours!
toplot <- toplot %>% mutate(r2 = purrr::pmap_dbl(list(c1=cpg1, c2=cpg2), check_cpg_r2))

# zoomed plot of r2s VS distances
zoomed <- toplot %>% ggplot(aes(x=distance, y=r2)) + geom_bin2d(bins=10000) + coord_cartesian(xlim=c(1,6e3))
# big picture plot
full <- toplot %>% ggplot(aes(x=distance, y=r2)) + geom_bin2d(bins=10000)
pg <- plot_grid(full, zoomed, ncol=2, labels=c("normal", "zoomed"))
save_plot("cpg_correlations_EPIC_full_and_zoomed.pdf", pg, ncol=2)

# create a plot with less bins
full2 <- toplot %>% ggplot(aes(x=distance, y=r2)) + geom_bin2d(bins=100)
save_plot("cpg_correlations_EPIC_full_bins100.pdf", full2, ncol=1)

write_tsv(path="results/current/EPIC_r2s_with_distances.tsv", toplot)




# ------------------------------------------------------------------------------
# Here we do a novel analysis, simply checking the fraction of all our sentinels
# with a proxy in EPIC

rm(list=ls())
gc(full=T)

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(reshape2)
library(parallel)
theme_set(theme_cowplot())
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
meth <- readRDS("results/current/methylation_residualized.rds")
sentinels <- read_tsv("data/benni.lehne/all_r02_110417_combined.tsv")
cpg_annot <- read_tsv("results/current/cpg_annotation_both_arrays.tsv")

# load cpg ranges for distance matching
cpg_ranges <- get_850k_cpg_ranges()

# get unique sentinel cpg ranges
sentinel_cpgs <- sentinels %>% pull(cpg.sentinel) %>% unique
sentinel_cpg_ranges <- cpg_ranges[sentinel_cpgs[sentinel_cpgs %in% names(cpg_ranges) &
                                                sentinel_cpgs %in% colnames(meth)]]

# Number of sentinels available on array:
print(length(sentinel_cpg_ranges))

# get epic cpg ranges
epic_cpgs <- filter(cpg_annot, is_epic & !is_450k) %>% pull(id)
epic_cpg_ranges <- cpg_ranges[epic_cpgs]

# ------------------------------------------------------------------------------
# for each sentinel cpg, get the EPIC CPG nearest to it
# ------------------------------------------------------------------------------
nearest_to_sentinel <- epic_cpg_ranges[nearest(sentinel_cpg_ranges,
                                               epic_cpg_ranges)]

# ------------------------------------------------------------------------------
# get all r^2s for the nearest entities
# ------------------------------------------------------------------------------
check_cpg_r2 <- function(c1, c2) {
  meth_r2 <- cor(meth[,c1], meth[,c2],
                                 use="complete")^2
  return(meth_r2)
}

threads <- 6

r2s <- unlist(mclapply(1:length(sentinel_cpg_ranges), function(i) {
  check_cpg_r2(names(sentinel_cpg_ranges)[i], names(nearest_to_sentinel)[i])
}, mc.cores = threads))

table(r2s>0.2)
#
#
#
