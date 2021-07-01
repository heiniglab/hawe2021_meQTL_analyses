#!/bin/bash
chr=$SLURM_ARRAY_TASK_ID
python3 -m tensorqtl data/current/kora/AffyAxiom/plink_methylation/KORAS4F4_methylation_chr${chr} results/current/tensorqtl/methylation_${chr}.bed.gz results/current/tensorqtl/interactions_bmi_chr${chr} --pval_threshold 1e-5 --covariates results/current/tensorqtl/covariates.txt --mode cis_nominal --interaction results/current/tensorqtl/covariates_utbmi.txt