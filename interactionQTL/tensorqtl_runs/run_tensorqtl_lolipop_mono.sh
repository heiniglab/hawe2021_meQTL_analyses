#!/bin/bash
chr=$SLURM_ARRAY_TASK_ID
python3 -m tensorqtl data/current/lolipop/plink/OmniEE/genotypes_chr${chr} results/current/tensorqtl/OmniEE/methylation_${chr}.bed.gz results/current/tensorqtl/OmniEE/interactions_Mono_chr${chr}  --maf_threshold_interaction 0.05 --pval_threshold 1e-5 --covariates results/current/tensorqtl/OmniEE/covariates.txt  --mode cis_nominal --interaction results/current/tensorqtl/OmniEE/covariates_Mono.txt
