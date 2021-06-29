## General
This project describes an meQTL pipeline for the KORA S4F4 genotypes and 850k methylation array.

The complete pipeline is implemented via Snakemake (version 5.8.2) and can be used to generate all meQTLs from the available raw data on the AME compute cluster (LISA).
The file `Snakefile` in the project base directory is the main workflow file containing and connecting all rules (i.e. directives/job definitions) to run the pipeline.

## Conda environment
It is recommended to run the pipeline in a pre-specified conda environment. As it is necessary to install the GitHub version of
matrix-EQTL, this step has to be done manually after creating the environment. In addition, the 850k array annotation package `IlluminaHumanMethylationEPICanno.ilm10b2.hg19` failed to install when using conda, so we have to install it manually, too.

First, install the conda environment found under `envs/rbio.yaml` (conda has to be installed on the system):

```
conda env create -f envs/rbio.yaml
``` 

> Note that this is not a 'minimal' environment, but rather a standard environment containing several 'nice to have' packages.

Once this is finished, activate the environment using

```
conda activate rbio
``` 

and start an R terminal.

Now we can install 1) devtools + MatrixEQTL and 2) `IlluminaHumanMethylationEPICanno.ilm10b2.hg19` (this will take a long while!):

```
# matrixEQTL
install.packages("devtools")
devtools::install_github("andreyshabalin/MatrixEQTL")

# annotation package
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
```

## Configuration
Configuration file can be found under `configs/workflow.json`. At the moment only the following options exist:

* `maf` : the minor allele frequency threshold to be applied when converting impute2 genotypes to dosages.
* `post_maf` : the minor allele frequency threshold to be applied after matrixEQTL has been run. useful to define different summaries.
* `cisDist` : the maximal distance between SNP-CpGs pairs until which it counts as a 'cis-pairs'. cis-pairs will be reported in a separate file.
* `pvThresholdCis` : The p-value threshold to be applied for cis-meQTLs. Only associations with `p-value < pvThresholdCis` will be reported
* `pvThresholdTrans` : The p-value threshold to be applied for trans-meQTLs. Only associations with `p-value < pvThresholdTrans` will be reported
* `significant_pv_cutoff` : For final results files: at which threshold should an association be called significant? e.g.: 1e-14

## Execution
To run the full pipeline, simply call snakemake with the default target (i.e. execute `snakemake` in the project directory).
In order to take advantage of multiple cores and get some additional outputs consider calling

```
./scripts/run_pipeline.sh <num-parallel-tasks>
```

> Note that at the moment you need to manually call the `all_dosage_to_matrixEQTL` and the `all_methylation_to_matrixEQTL` rules **before** running the line above!

The above script is just a simple wrapper for running snakemake with parallel task processing enabled and also restricted to a maximal amount of 160GB RAM.
It will further generate three distinct visualizations (as graphs)  describing the current workflow.

It's tuned such that it works well with the maximum available VM on LISA (80 cores, 160GB mem), and the script `[scripts/run_pipeline.sh](scripts/run_pipeline.sh)` (see above)
is available to get it started easily for this configuration.

Expected runtime for running matrixEQTL on all inputs is about 2 days, for this it will run 6,396 matrixEQTL jobs, each specified with 12,000 Gbyte Ram and 6 threads.
A single one of those jobs takes about 3 minutes.

## Covariates 

We used the following model for obtaining EPIC meQTLS in the first rnu:

```
meth ~ geno + sex + age + bmi + wbc + PCs1..20 + houseman estimates (except gran) + batch
```

In a second run, in order to be closer to the original model definition of the paper,
we used the following covariates:

```
meth ~ geno + sex + age + bmi + wbc + houseman estimates (except gran) + plate + batch
```

