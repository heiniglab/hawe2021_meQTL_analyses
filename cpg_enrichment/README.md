# Functional characterization of meQTL CpG sites
This file includes description for the following R programs, (1) univariate association analysis which tests associations of individual CpG sites with other traits and (2) permutation analysis, which tests for enrichment of associations in a select set of CpG sites (meQTL CpGs) vs a background set of CpG sites.

## Univariate association analysis code

  The code in [Univariate_Code.R](Univariate_Code.R) performs univariate association analysis between a single trait (or phenotype) and all CpG markers (provided by the user). 
  Function "lm()" was used for quantitative traits. Function "glm()" was used for categorical traits. 
  It can be implemented in parallel (refer below for details).
  It doesn't require installation of any additional R packages

### Input files and their structure

  For association testing, this code requires two input files, 
  1. Phenotype (or traits) data ( e.g. bmi, T2D status, metabolite concentrations etc. ). 
      In the current implementation, we have provided this data as a binary RData file. 
      Each row in this file represents data for 1 sample and each column represents a phenotype or life style trait.
      Covariates used in the regression model are also included in this data file (e.g. age, gender, cell type composition etc. ) 
      At least one column should contain a sample identifier, which should also be present in the CpG data file, to correctly 
      match the CpG data and phenotype data between samples. In our data set, sample ID column header was "X450K_Array_ID"
  2. CpG marker data
      Each row represents methylation levels for a single CpG marker across samples. 
      Each column represent methylation levels of all CpG marker for a single sample.
      CpG data was call rate filtered and quantile normalized prior to this analysis (please refer to method description in the manuscript)

### Parallel implementation

  Phenotype data file can contain data for multiple phenotypes (i.e.in different columns) 
  This code can be implemented in parallel in compute clusters. Each node will perform the analysis for 1 phenotype column 
  For parallel implementation, user should provide a single input argument referring to the phenotype column index in the data file    


### Steps included in the code

  1. loading input files
  2. extracting the correct batch of samples, batch 1 in this code.  This will depend on the input data
  3. selecting the correct phenotype column in the data based on user provided input argument
  4. remove samples with missing (NA) values in phenotype data 
  5. filter out samples missing in any one of the phenotype file or CpG input data files  
  6. log transform metabolite concentration values
  7. scale CpG values
  8. subset large CpG data into small subsets and perform regression analysis
  9. write regression analysis results to output files

### Output files

  Output files are named using phenotype names and are in .csv format
  Output contains univariate association results from the regression analysis 


## Permutation analysis code

 The code in [Permutation_Testing.R](Permutation_Testing.R) generates randomly sampled sets of Sentinel CpG markers and matching Background CpG markers. 
 For a user specified univariate association results data ( single phenotype vs. CpG levels ), this code performs permutation analysis.
 permutation analysis can be performed for any one of the three categoires of CpGs (a) cis, (b) long-range and (c) trans
 User can adjust (a) number of randomly sampled sets and (b) number of CpGs per sampled set. By default, 1000 sets are sampled, each with 1000 CpGs in it.

### Input files and user specified parameters

 1. category tag: 'cis' or 'long-range' or 'trans'
 2. phenotype name 
 3. univariate results data file for the chosen phenotype
 4. file containing the full list of Sentinel CpGs. This file should have a column specifying whether the CpG is 'cis' or 'long-range' or 'trans'
 5. file containing the full list of cosmopolitan CpGs
 6. file containing mean and SD values for each CpG marker


### Steps included in the code

 1. read the full Sentinel CpG list and extract separately the correct category of CpGs (i.e. 'cis' or 'long-range' or 'trans')
 2. read list of cosmopolitan CpGs and generate Background CpGs set from Cosmopolitan set
 3. generate randomly sampled sets of Sentinel and Baackground CpGs
 4. perform permutation analysis
 5. write permutation analysis results to output files. Save randomly sampled CpG sets to separate files

### Output files

 1. generated random samples of Sentinel and Background CpGs are saved as text files
 2. permutation analysis results are saved as .csv files

