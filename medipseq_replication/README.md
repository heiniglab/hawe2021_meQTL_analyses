# medipseq_replication

## Background

This code runs analysis for the paper "Genetic variation influencing DNA methylation provides new insights into the molecular pathways regulating genomic function" by Hawe et al. Specifically, it runs the analysis for replicating the HSM-LD block associations from the MeDIP-seq results of Bell et al. (PMID: 29295990).

## Code

The R scripts are found in 'code/' and are the two files: 

* *step_1_map_HSM_LD_assoc_to_SNP_CpG.r*

We take the HSM-LD block associations from PMID: 29295990 and map the associations to pairs of SNPs and CpGs analyzed in our data. 

* *step_2_calculate_enrichment_above_bg.R*

Based on summary statistics, we see how many of the HSM-LD block associations we can replicate with pairs of SNPs and CpGs in our data.  We check also the summary statistics for the 100 sets of matched SNP-CpG pairs (background) to see if the HSM-LD-mapped SNP-CpG pairs show a greater number of significant hits than the background. 

## Data

The data required for the code are found in 'data/' and are the following:

* *mergeFinal_hsm_LDblock_GWASSNPs.txt*

Results from PMID: 29295990.

* *illumina_450K_annot_pruned.txt*

Pruned Illumina annotation file. 

* *kora_cpgs.RData*

CpG sites used in our analysis in KORA (ie, passing QC).

* *kora_snps.RData*

SNPs used in our analysis in KORA (ie, passing QC).

* *medip_replication_results.txt*

The linear regression (association) results of the HSM-LD-mapped SNP-CpG pairs.

* *medip_matched_background_results.txt*

The linear regression (association) results of the 100 sets of matched SNP-CpG pairs (background).

## Notes

It is possible to completely reconstruct the results from the data found here. However, the data presented here are not raw data. Specifically, after step_1, the following cannot be performed due to data privacy reasons: 

* we mapped the CpG-SNP pairs generated in step_1 to 100 sets of matched background CpG-SNP pairs
* for the the CpG-SNPs pairs and then their matched sets, we used individual-level KORA F4 data to check the associations, using linear regression

The results of these analyses (*medip_replication_results.txt* and *medip_matched_background_results.txt*) are then used as inputs for step_2.
