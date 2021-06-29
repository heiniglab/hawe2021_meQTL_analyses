#' -----------------------------------------------------------------------------
#' Check MAF distribution for cis/longrange/trans meQTLs
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jan  07 14:32:38 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries, set parameters.")
# ------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
# inputs
fannot <- snakemake@input$annot
fcis <- snakemake@input$cis
ftrans <- snakemake@input$trans
flongrange <- snakemake@input$longrange

# outputs
fout_plot <- snakemake@output$plot
fout_cis <- snakemake@output$cis
fout_longrange <- snakemake@output$longrange
fout_trans <- snakemake@output$trans
fout_combined <- snakemake@output$combined

# params
maf <- snakemake@params$maf
print(paste0("MAF is set to: ", maf, "."))

# ------------------------------------------------------------------------------
# join MAF info to meQTL results
# ------------------------------------------------------------------------------
annot <- readRDS(fannot)
cis <- readRDS(fcis)
longrange <- readRDS(flongrange)
trans <- readRDS(ftrans)

cis_maf <- inner_join(cis, annot, by = c("SNP" = "rsid"))
trans_maf <- inner_join(trans, annot, by = c("SNP" = "rsid"))
longrange_maf <- inner_join(longrange, annot, by = c("SNP" = "rsid"))

# ------------------------------------------------------------------------------
# check total number of associations for adjusted MAF (i.e. MAF > 0.5 -> MAF = 1-MAF)
# ------------------------------------------------------------------------------
filter_maf <- function(data, title, plot = T, maf=0.01) {
  print(paste0("Checking ", title))  
  print(paste0("Number of associations before filtering: ", nrow(data)))
  
  adjusted <- data %>% mutate(MAF = case_when(MAF > 0.5 ~ 1-MAF, 
                                              MAF < 0.5 ~ MAF)) %>%
  filter(MAF > maf)
  
  print(paste0("Number of associations after filtering: ", nrow(adjusted)))
  
  if(plot) {
    # plot simple histogram
    adjusted %>% 
      ggplot(aes(x=(MAF))) + 
      geom_histogram(bins=100) + 
      background_grid(major="xy") + 
      labs(title=title)
  } else {
    adjusted
  }
}

gp_cis <- filter_maf(cis_maf, "cis meQTLs", maf=maf)
gp_longrange <- filter_maf(longrange_maf, "longrange meQTLs", maf=maf)
gp_trans <- filter_maf(trans_maf, "trans meQTLs", maf=maf)

# plot everything in single grid
save_plot(file_name=fout_plot, 
          plot=plot_grid(gp_cis, gp_longrange, gp_trans, ncol=3, nrow=1), 
          ncol=3)

#' -----------------------------------------------------------------------------
print("Saving refiltered associations.")
#' -----------------------------------------------------------------------------
new_cis <- filter_maf(cis_maf, "cis meQTLs", F, maf)
new_longrange <- filter_maf(longrange_maf, "longrange meQTLs", F, maf)
new_trans <- filter_maf(trans_maf, "trans meQTLs", F, maf)

saveRDS(new_cis, fout_cis)
saveRDS(new_longrange, fout_longrange)
saveRDS(new_trans, fout_trans)

combined <- bind_rows(new_cis, new_longrange, new_trans)
saveRDS(combined, fout_combined)

#' -----------------------------------------------------------------------------
print("SessionInfo:")
#' -----------------------------------------------------------------------------
sessionInfo()
