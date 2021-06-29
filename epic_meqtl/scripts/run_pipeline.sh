#!/bin/bash

## execute line below BEFORE executing the script to
## make snakemake and all beeded packages available
#conda activate rbio

## number of parallle tasks in snakemake (e.g. available cores..)
par=$1
if [[ -z $par ]] ; then par=80 ; fi
 
## call snakemake using the main target rule
## specify a maximum of ${par} cores and 160GB RAM (current limits of VMs)
# export openblas thread definition (R takes #CORES as default..)
export OPENBLAS_NUM_THREADS=1 
#nohup snakemake -k --resources mem_mb=163000 -j ${par} all &> all.out
nohup snakemake --nolock -k --resources mem_mb=163000 -j ${par} merge_significant &> merge_significant.out && \
  ## create overview graphs for workflow
  snakemake --dag results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/matrixEQTL_outputs/chr20_meth1_cis_associations.out | dot -Tpdf > dag_chr20.pdf && \
  snakemake --filegraph results/current/maf0.01_cis1_trans1e-5_sig5e-2_postMAF0.01/matrixEQTL_outputs/chr20_meth1_cis_associations.out | dot -Tpdf > filegraph_chr20.pdf && \
  snakemake --rulegraph all | dot -Tpdf > rulegraph.pdf 
