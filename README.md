# Code for the meQTL analyses presented in Hawe et al. 2021 Nature Genetics

Citation:
```
Hawe et al 2021 Nature Genetics
```

Table of Contents
=================

   * [Identification of meQTLs](#identification-of-meqtls)
   * [Identification of eQTMs](#identification-of-eqtms)
   * [Functional analyses of meQTLs](#functional-analyses-of-meqtls)
      * [Chromatin state enrichment](#chromatin-state-enrichment)
      * [HiC enrichment](#hic-enrichment)
      * [Regulator enrichment](#regulator-enrichment)
      * [Enrichment of transcription factor binding sites at trans-meQTL CpG sites](#enrichment-of-transcription-factor-binding-sites-at-trans-meqtl-cpg-sites)
      * [Network analysis](#network-analysis)
   * [Identification of iQTL](#identification-of-iqtl)
   * [EPIC comparison](#epic-comparison)

## Identification of meQTLs

## Identification of eQTMs

## Functional analyses of meQTLs

### Chromatin state enrichment

### HiC enrichment

All details for the HiC contact enrichment are provided in a separate [README](hic_enrichment/README.md).

### Regulator enrichment

Enrichment of (epigenetic) regulator sets provided in the Lemire et al. 2015 Nat Comms paper, see [regulator_enrichment/](./regulator_enrichment/)

### Enrichment of transcription factor binding sites at trans-meQTL CpG sites

We test for over-representation of TFBS of specific factors at the CpG sites that share trans meQTL. To reproduce the results in Figure 3 and the extended data figure 4 follow these steps.

First run `R/annotate-cpgs.R` this creates CpG contexts of sizes 2, 100, 500,
1000, 5000, 10000

```{bash}
for i in 2 100 500 1000 5000 10000; do
  qsub -V -cwd -l job_mem=15G -N tfbs_context_$i -q long@@core R/tfbs-enrichment.R --size $i --resample 0
done
```

Also run the TFBS analysis on the union of binding sites for each factor over
the different experiments

```{bash}
qsub -V -cwd -l job_mem=15G -N tfbs_context_$i -q long@@core R/tfbs-enrichment.R --size 100 --resample 0 --prefix results/current/enrichment_chipseq_context_100_resample_0_by_tf --tfbs results/current/cpgs_with_chipseq_context_100_by_tf.RData
```

To run the analysis with resampling instead of the analytical P-values you can specify the `--resample<int>` option and specify how many resampling runs should be done to determine empirical P-values. 

### Network analysis

Code for the network analysis based on random walks is provided in the R package [QTLnetwork](https://github.com/heiniglab/QTLnetwork).

## Identification of iQTL

## EPIC comparison

We investigated to which degree the increased coverage and resolution of the EPIC array could impact our conclusions from the 1) cis regulatory element enrichment and 2) TFBS enrichment analyses.
The scripts to calculate and compare the enrichments for the original 450k and the new EPIC meQTL pairs are deposited under [epic_comparison/](epic_comparison/)
