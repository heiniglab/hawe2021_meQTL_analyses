# Code for the meQTL analyses of Hawe et al. 2021 Nature Genetics

Citation:

**Genetic variation influencing DNA methylation provides new insights into the molecular pathways regulating genomic function.**

Johann S Hawe, Rory Wilson, Katharina Schmid, Li Zhou, Lakshmi Lakshmanan, Benjamin C Lehne, Brigitte Kühnel, William R Scott, Matthias Wielscher, Yik Weng Yew, Clemens Baumbach, Dominic P Lee, Eirini Marouli, Manon Bernard, Liliane Pfeiffer, Pamela Matías-García, Matias I Autio, Stephane Bourgeois, Christian Herder, Ville Karhunen, Thomas Meitinger, Holger Prokisch, Wolfgang Rathmann, Michael Roden, Sylvain Sebert, Jean Shin, Konstantin Strauch, Weihua Zhang, Wilson LW Tan, Stefanie M. Hauck, Juliane Merl-Pham, Harald Grallert, Eudes GV Barbosa, MuTHER Consortium, Thomas Illig, Annette Peters, Tomas Paus, Zdenka Pausova, Panos Deloukas, Roger SY Foo, Marjo-Riitta Jarvelin, Jaspal S Kooner, Marie Loh, Matthias Heinig, Christian Gieger, Melanie Waldenberger, and John C Chambers.


## Table of Contents


   * [Identification of meQTLs](#identification-of-meqtls)
      * [EPIC meQTL](#epic-meqtl)
   * [meQTL replication meDIP-seq](#meqtl-replication:-medip-seq) 
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

### EPIC meQTL

We ran a smaller analysis to evaluate the potential impact of the higher coverage and density of the EPIC array on our analyses.
We ran a separate meQTL pipeline using the available KORA EPIC data, details can be found under [epic_meqtl/](epic_meqtl/).

## meQTL replication: meDIP-seq

Identified meQTL were corroborated using several replication analyses.
The meDIP-seq replcation is provided in a separate repository: [meDIP-seq replication](https://github.com/Rory-R-Wilson/medipseq_replication/)

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

We tested if meQTL relationships might be influenced by environmental context. For this, we ran interaction analyses for SNP-CpG pairs using linear regression models with an interaction term between SNP and phenotype. 

```
meth ~ geno + pheno + geno:pheno + covariates
```

The following phenotypes were examined: smoking (yes/no), BMI and estimated proportions of CD8 T cells, CD4 T cells and monocytes. The analysis was once performed restricted to the cosmopolitan meQTLs and once in a genome-wide cis analysis.

### iQTLs at cosmopolitan meQTLs

### Genome-wide cis-iQTLs

For the genome-wide cis iQTL analysis, we ran the analysis using [tensorqtl](https://github.com/broadinstitute/tensorqtl). We performed the analysis separetly for each cohort and phenotype and parallelized over the chromosomes, submitting to a GPU server with a slurm workload manager (see scripts in `interactionQTL/tensorqtl_runs`), e.g.

```{bash}
sbatch --array=1-22 interactionQTL/tensorqtl_runs/run_tensorqtl_kora_bmi.sh 
```
Afterwards, significant iQTLs were determined using a Bonferroni corrected multiple testing threshold of 0.05 with

```{bash}
python -u interactionQTL/tensorqtl_runs/filter_iQTLs.py
```

## EPIC comparison

We investigated to which degree the increased coverage and resolution of the EPIC array could impact our conclusions from the 1) cis regulatory element enrichment and 2) TFBS enrichment analyses.
The scripts to calculate and compare the enrichments for the original 450k and the new EPIC meQTL pairs are deposited under [epic_comparison/](epic_comparison/)
