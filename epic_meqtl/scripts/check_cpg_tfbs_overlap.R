library(rtracklayer)
library(data.table)
source("scripts/lib.R")


# get tfbs annotation for relevant CpGs
tfbs.ann <- load_tfbs_annot(fread("results/current/cpgs_to_check_with_tfbs_50.txt", header=F, data.table=F)$V1, 100)

# how many cpgs have any TFBS?
ncpgs_with_tfbs <- sum(apply(tfbs.ann, 1, any))
ncpgs_with_tfbs

cpgs_with_tfbs <- rownames(tfbs.ann[apply(tfbs.ann, 1, any),])
write_tsv(tibble(cpgs=cpgs_with_tfbs), col_names=F, path="cpgs_with_tfbs.txt")

# get fraction of r2>0.2 (with nearest sentinel) CpGs with trans-meQTL:
# bash: 
# awk '{if($5+0 < 1e-14) print;}' results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/meqtls/trans_associations_chr* | \
  fgrep -w -f results/current/cpgs_to_check_with_tfbs.txt | cut -f 2 | sort | uniq > trans_cpgs_to_check.txt
or (for CpGs with TFBS):
  fgrep -w -f cpgs_with_tfbs.txt | cut -f 2 | sort | uniq > trans_cpgs_with_tfbs.txt

