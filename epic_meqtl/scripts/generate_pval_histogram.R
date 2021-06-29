#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

# we only expect the matrix EQTL output file
fmeqtls <- commandArgs(trailingOnly = T)[1]

# load data
meqtls <- read_tsv(fmeqtls)

# create histogram plot
meqtls %>% ggplot(aes(x=`p-value`)) + 
  geom_histogram(bins=100) +
  background_grid(major="xy")
