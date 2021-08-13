## README

Enrichment analyses for the lists of epigenetic regulators provided in the Lemire et al. 2015 paper (https://www.nature.com/articles/ncomms7326)

Lists of regulators are provided in lemire2015/ sub-directory.

Other input data are provided via [zenodo](https://zenodo.org/record/5196216#.YRZ3TfJxeUk)

Can be run via the provided snakemake file, e.g.:

```{bash}
snakemake --profile profiles/default regulator_enrichment
```
