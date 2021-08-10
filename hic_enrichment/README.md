## README

Input files for this analysis can be obtained from [zenodo](ZENODO_LINK).

The `Snakefile` in this directory contains all rules used to evaluate the HiC
enrichment for `longrange` (TAD + HiC) and `trans` meQTL pairs.
Can be executed using [snakemake](https://snakemake.readthedocs.io), e.g.:

```{bash}
snakemake --profiles profiles/default all_hic
```

The `report.Rmd` can be generated after running the snakemake workflow.
The markdown was originally in a separate 'master-report' for the full meQTL
project, and the parts relevant for the HiC enrichment have been extracted.
Generate the report e.g. using:

```{bash}
R -e "rmarkdown::render('report.Rmd')"
```