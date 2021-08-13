# Code for the meQTL analyses of Hawe et al. 2021 Nature Genetics

Citation:

**Genetic variation influencing DNA methylation provides new insights into the molecular pathways regulating genomic function.**

Johann S Hawe, Rory Wilson, Katharina Schmid, Li Zhou, Lakshmi Lakshmanan, Benjamin C Lehne, Brigitte Kühnel, William R Scott, Matthias Wielscher, Yik Weng Yew, Clemens Baumbach, Dominic P Lee, Eirini Marouli, Manon Bernard, Liliane Pfeiffer, Pamela Matías-García, Matias I Autio, Stephane Bourgeois, Christian Herder, Ville Karhunen, Thomas Meitinger, Holger Prokisch, Wolfgang Rathmann, Michael Roden, Sylvain Sebert, Jean Shin, Konstantin Strauch, Weihua Zhang, Wilson LW Tan, Stefanie M. Hauck, Juliane Merl-Pham, Harald Grallert, Eudes GV Barbosa, MuTHER Consortium, Thomas Illig, Annette Peters, Tomas Paus, Zdenka Pausova, Panos Deloukas, Roger SY Foo, Marjo-Riitta Jarvelin, Jaspal S Kooner, Marie Loh, Matthias Heinig, Christian Gieger, Melanie Waldenberger, and John C Chambers.


## Table of Contents


   * [Identification of meQTLs](#identification-of-meqtls)
      * [EPIC meQTL](#epic-meqtl)
      * [Conditional analysis and pruning](#conditional-analysis-and-pruning)
   * [meQTL replication meDIPseq](#meqtl-replication-medipseq) 
   * [Identification of eQTLs and eQTMs](#identification-of-eqtls-and-eqtms)
   * [Functional analyses of meQTL CpGs](#functional-analyses-of-meqtl-cpgs)
   * [Functional analyses of meQTL pairs](#functional-analyses-of-meqtl-pairs)
      * [Chromatin state enrichment](#chromatin-state-enrichment)
      * [HiC enrichment](#hic-enrichment)
      * [Regulator enrichment](#regulator-enrichment)
      * [Enrichment of transcription factor binding sites at trans-meQTL CpG sites](#enrichment-of-transcription-factor-binding-sites-at-trans-meqtl-cpg-sites)
      * [Network analysis](#network-analysis)
   * [Identification of iQTL](#identification-of-iqtl)
   * [EPIC comparison](#epic-comparison)
   * [Docker](#docker)

## Identification of meQTLs

### EPIC meQTL

We ran a smaller analysis to evaluate the potential impact of the higher coverage and density of the EPIC array on our analyses.
We ran a separate meQTL pipeline using the available KORA EPIC data, details can be found under [epic_meqtl/](epic_meqtl/).

### Conditional analysis and pruning
Conditional meQTL analysis and LD pruning is described here: [conditional_analysis/](conditional_analysis/).

## meQTL replication: meDIPseq

Identified meQTL were corroborated using several replication analyses.
The meDIP-seq replcation is provided in a separate repository: [meDIP-seq replication](medipseq_replication/)

## Identification of eQTLs and eQTMs

We computed both associations between genotypes and expression (expression ~ genotypes), commonly referred to as expression quantatiative trait loci (eQTL) and associations between between methylation and expression (expression ~ methylation), commonly referred to as expression quantiative trait methylation (eQTM).

The eQTL analysis is restricted to sentinel meQTL SNPs. Following the snakemake workflow [eqtl_sentinel_snp/Snakefile](eqtl_sentinel_snp/Snakefile.sm), eQTLs are calculated separately for each cohort (adjusting the expression for covariates and houseman groups before) and then the results are combined in a meta-analysis. 

A similar workflow is applied for the eQTMs in [eqtm_identification/Snakefile](eqtm_identificaton/Snakemake.sm), with a separate analysis per cohort followed by a meta-analsis. The eQTM analysis is carried out with two different versions of the expression and methylation data, once adjusted for only the genotype (files called GTregressed) and once corrected for both the genotype and cell proportions (files called CPregressed). The genotype adjustment is done using cis eQTL SNPs for the expression adjustment and using cis meQTL SNPs for the methylation adjustment. Additionally, the replication rates for the eQTMs between cohorts are calculated, both for Bonferroni and FDR corrected significance thresholds. The FDR calculation requires sorting the files locally with sufficient amount of RAM, this part is not incorporated in the snakemake workflow (please take a look at the documentation of the associated rule "replicate_eqtm" for details).

The analysis requires genotype, expression, methlyation and covariate data (including Houseman estimates of the cell type proportions) of both cohorts as input. Required R packages for the snakemake workflows are matrixeQTL, RhpcBLASctl, tidyverse, data.table, GenomicRanges, purrr (all included  in the [conda environment](epic_meqtl/envs/rbio.yaml)).

## Functional analyses of meQTL CpGs
CpG sites with meQTL were tested for association with other traits in our cohort. Based on these results meQTL CpG sites were compared to random background CpG sets to test for enrichment of trait associations. Analyses are described here: [cpg_enrichment/](cpg_enrichment/).

## Functional analyses of meQTL pairs

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

The following phenotypes were examined: smoking (yes/no), BMI and estimated proportions of CD8 T cells, CD4 T cells and monocytes. 

For the genome-wide cis iQTL analysis, we ran the analysis using [tensorqtl](https://github.com/broadinstitute/tensorqtl). We performed the analysis separetly for each cohort and phenotype and parallelized over the chromosomes, submitting to a GPU server with a slurm workload manager (see scripts in `interactionQTL/tensorqtl_runs`), e.g.

```{bash}
sbatch --array=1-22 interactionQTL/tensorqtl_runs/run_tensorqtl_kora_bmi.sh 
```
Afterwards, significant iQTLs were determined using a Bonferroni corrected multiple testing threshold of 0.05 with

```{bash}
python -u interactionQTL/tensorqtl_runs/filter_iQTLs.py
```
The analysis requires genotype, methlyation and covariate data of both cohorts as input.

## EPIC comparison

We investigated to which degree the increased coverage and resolution of the EPIC array could impact our conclusions from the 1) cis regulatory element enrichment and 2) TFBS enrichment analyses.
The scripts to calculate and compare the enrichments for the original 450k and the new EPIC meQTL pairs are deposited under [epic_comparison/](epic_comparison/)

## Docker

We provide some additional information as well as our [DOCKERFILE](docker/DOCKERFILE) under [docker/](docker/) in this repository.
