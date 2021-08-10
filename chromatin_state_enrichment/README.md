## cis-meQTL chromHMM enrichment

We investigated whether cis-meQTLs are enriched for enhancer/promoter chromHMM states
as compared to a random, matched background (150 iterations).

The main R-script for this analysis is located under [snp_cpg_annotation_overlap_batchjobs.R](../R/cre_enrichment/snp_cpg_annotation_overlap_batchjobs.R).

Input files for this analysis can be obtained from [zenodo](ZENODO_LINK).

In order to calculate all iterations, a HPC should be used since these are
some intensive calculations. In our case, we use the *all_snp_cpg_annotation_overlap*
target rule and make a cluster call like so:

```{bash}
nohup nice snakemake --profile profiles/slurm all_snp_cpg_annotation_overlap &
```

> NOTE: you need to create the 'slurm' snakemake profile for your environment for the 
above command to execute successfully.

The results are then finalized in the [report.Rmd](report.Rmd) markdown under "cis meQTL chromHMM enrichment".
The report also contains a brief evaluation/comparison of the EPIC based meQTL results.
