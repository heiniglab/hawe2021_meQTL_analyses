## README

Input files for this analysis can be obtained from [zenodo](https://zenodo.org/record/5196216#.YRZ3TfJxeUk).

### TFBS enrichment TODO
TFBS enrichment via Snakemake, see [Snakefile](Snakefile).

### CRE enrichment TODO

Cis-regulatory element enrichment uses these scripts:

* Overlap in EPIC pairs: `R/cre_enrichment/annotation_overlap_EPIC.R`
* Overap for background pairs: `R/cre_enrichment/annotation_overlap_EPIC_background.R`
* Helper methods for enrichment: `R/cre_enrichment/annotation_overlap_methods.R`
