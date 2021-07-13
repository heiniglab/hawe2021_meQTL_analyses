# Code description: 
# ------------------
#  This code performs univariate association analysis between a single trait (or phenotype) and all CpG markers (provided by the user). 
#  Function "lm()" was used for quantitative traits. Function "glm()" was used for categorical traits. 
#  It can be implemented in parallel (refer below for details).
#  It doesn't require installation of any additional R packages

# Input files and their structure
# --------------------------------
#  For association testing, this code requires two input files, 
#  (1) Phenotype (or traits) data ( e.g. bmi, T2D status, metabolite concentrations etc. ). 
#      In the current implementation, we have provided this data as a binary RData file. 
#      Each row in this file represents data for 1 sample and each column represents a phenotype or life style trait.
#      Covariates used in the regression model are also included in this data file (e.g. age, gender, cell type composition etc. ) 
#      At least one column should contain a sample identifier, which should also be present in the CpG data file, to correctly 
#      match the CpG data and phenotype data between samples. In our data set, sample ID column header was "X450K_Array_ID"
#  (2) CpG marker data
#      Each row represents methylation levels for a single CpG marker across samples. 
#      Each column represent methylation levels of all CpG marker for a single sample.
#      CpG data was call rate filtered and quantile normalized prior to this analysis (please refer to method description in the manuscript)

# Parallel implementation
# ------------------------
#  Phenotype data file can contain data for multiple phenotypes (i.e.in different columns) 
#  This code can be implemented in parallel in compute clusters. Each node will perform the analysis for 1 phenotype column 
#  For parallel implementation, user should provide a single input argument referring to the phenotype column index in the data file    


# Steps included in the code
# --------------------------
#  (1) loading input files
#  (2) extracting the correct batch of samples, batch 1 in this code.  This will depend on the input data
#  (3) selecting the correct phenotype column in the data based on user provided input argument
#  (4) remove samples with missing (NA) values in phenotype data 
#  (5) filter out samples missing in any one of the phenotype file or CpG input data files  
#  (6) log transform metabolite concentration values
#  (7) scale CpG values
#  (8) subset large CpG data into small subsets and perform regression analysis
#  (9) write regression analysis results to output files


# Output files
# ------------
#  Output files are named using phenotype names and are in .csv format
#  Output contains univariate association results from the regression analysis


# code last updated: 20/Nov/2020
# comments added later, without modifying the code 


  cat('\f')
  
  rm(list = ls())
  
  graphics.off()
  
  args=(commandArgs(trailingOnly = TRUE))   # for parallel implementation, user need to specify 1 input argument
  
 
  
  # load input data 
  # -----------------
  
  
  load("/home/projects/12000713/Lakshmi/Lolipop_Phe_Metabolomics.RData") # Metabolomics_Phe_DF                   # provide correct file directories 
  
  phe = data.frame()
  
  phe = Metabolomics_Phe_DF
  
  rm(list=c('Metabolomics_Phe_DF'))
  
  Batch_Tag_List = c('1')    # batch 1 is the largest batch in our data set, analysis was performed for batch 1 samples  
  
  col_number_val = as.integer( args[1] )   # phenotype column specified by user in input argument
  
  # in our dataset, metabolic phenotypes start from column 133
  # to accomodate for this,  if user specifies metabolite 1,  then it will be column (132 + 1) in our input data
  
  metabolite_col_number = 0;            metabolite_col_number = col_number_val + 132   
  
  Header_List = names(phe)
  
  curr_metab_name = '';                 curr_metab_name = Header_List[metabolite_col_number]
  
  for (batch_tag in Batch_Tag_List) {   # for loop to analyze samples from multiple batches, but here we use only batch 1 
    
    load( paste("/home/projects/12000713/Lakshmi/Sep23_2020_QN/Batch", batch_tag, "_QN_CG_Data.RData", sep="") )   # provide correct file directories
    
    IA_epimigrant_beta = CG_Data_QN
  
    rm(list = c('CG_Data_QN'))
  
    EpiData_Header_Names = colnames(IA_epimigrant_beta)
  
  
    # Batch 
    # ---------------

    Batch1_Indices = c();               Batch1_Indices = which( (phe$Meth_Batch == as.numeric(batch_tag) ) )   # extracting batch 1 samples
  
    Batch1_Phe = data.frame();          Batch1_Phe = phe[Batch1_Indices,]
  
  
    # # in our data, zero values represent NA values (i.e. missing values)
    # # separate preliminary analysis was done for each phenotype, to assess their data values distribution 
    # # Here we identify and remove samples with zero values for phenotypes 
    # # -------------------------------------------------------
                                                                                                               
    Zero_Indices = c();                 Zero_Indices = which(Batch1_Phe[,metabolite_col_number] == 0)          
                                                                                                               
    if( length(Zero_Indices) > 0 ) {
  
      Batch1_Phe = Batch1_Phe[-c(Zero_Indices), ]
  
    }
  
    # # --------------------------------------------------------
  
    Batch1_ArrayIDs = as.character( Batch1_Phe$X450K_Array_ID )
  
    # for performing regression analysis we need both phenotype data and CpG data
    # since phenotype data and CpG data are obtained as separate files, there could be some samples that do not contain both phenotype and CpG data
    # hence, we check for sample IDs, which are present in phenotype data but not present in CpG data
    
  
    CR_Filtered_Samples = setdiff(Batch1_ArrayIDs, EpiData_Header_Names)
  
    Filtered_Sample_Indices = c()
  
    for(filter_loop in 1:length(CR_Filtered_Samples)){
  
      curr_index = 0;                   curr_index = which(Batch1_ArrayIDs == CR_Filtered_Samples[filter_loop])
    
      Filtered_Sample_Indices = append(Filtered_Sample_Indices, curr_index)
    
    } # filter loop
  
    if( length(Filtered_Sample_Indices) > 0 ) {
    
      Batch1_Phe = Batch1_Phe[-c(Filtered_Sample_Indices), ]
    
    }
  
    Batch1_ArrayIDs = c()
  
    Batch1_ArrayIDs = as.character( Batch1_Phe$X450K_Array_ID )
  
  
    # in CpG data, each column represents a sample. CpG data contains samples from different batches.
    # here, we identify the column indices of batch 1 samples within the whole CpG data
    
  
    Batch1_Meth_Col_Indices = c()
  
    for ( sample_loop in 1:nrow(Batch1_Phe) ) { # 1:nrow(Batch1_Phe)
    
      sample_column_index = 0;         sample_column_index = which( EpiData_Header_Names == Batch1_ArrayIDs[sample_loop] )
    
      # #if(length(sample_column_index) == 1) {
    
      Batch1_Meth_Col_Indices = append( Batch1_Meth_Col_Indices, sample_column_index )
    
      # # }
    
    }
  
  
    # to avoid performing analysis on all ~475,000 CpGs together, we divide total CpGs into smaller subsets (50 subsets with 9,500 CpGs) 
    
    CG_unit_value = 9500
  
    max_CG_value = nrow(IA_epimigrant_beta)
  
    total_jobs = 50      
  
    t1 = proc.time()
  
   
    for( job_number in 1:total_jobs ){ # 1:total_jobs
    
      CG_start_row = 0;                   CG_start_row =  ( (job_number-1) * CG_unit_value ) + 1             # 1 # as.integer( args[1] )
    
      CG_stop_row = 0;                    CG_stop_row = ( (job_number-1) * CG_unit_value ) + CG_unit_value   # 1000 #as.integer( args[2] )
    
      if( job_number == 50 ) {
      
        CG_stop_row = max_CG_value
      
      }  
      
      # # result files are named using phenotype name (metabolite name in this specific case)
    
      results_filename = paste( curr_metab_name, '_CG_B', batch_tag, '_Results_Job', job_number, '.csv', sep='')         
                                                                                                                        
    
    
      # dependent variable, independent variables preparation
      # Columns and rows are somewhat reversed in the input Phenotype data and CpG data 
      # (e.g. in phenotype data each row is one sample, in CpG data each column is 1 sample)
      # CpG data is transposed, so that it will be easier to merge it with phenotype data
    
      Selected_CG_Data = data.frame();            Selected_CG_Data = IA_epimigrant_beta[CG_start_row:CG_stop_row,]
    
      total_markers = nrow(Selected_CG_Data)
    
      Selected_CG_Data = Selected_CG_Data[1:total_markers, Batch1_Meth_Col_Indices]
    
      #print( ncol(Selected_CG_Data) )
    
      Independent_Variables_Data =  as.data.frame( t( Selected_CG_Data ) )          
    
      Independent_Variables_List_CG = names(Independent_Variables_Data) # CG marker IDs
    
      Independent_Variables_Data$X450K_Array_ID = row.names(Independent_Variables_Data)
    
    
      # Reordering the methylation data through 'merge' function
    
      Updated_Phe_Data = merge(Batch1_Phe, Independent_Variables_Data, by.x = 'X450K_Array_ID', all = TRUE, sort = FALSE)
    
      Updated_Phe_Data$Sex = as.factor(Updated_Phe_Data$Sex)      # assigning Sex as a categorical variable 
    
      total_samples = 0;                          total_samples = nrow(Updated_Phe_Data)
    
    
      # # log transform metabolite concentrations
      # # ----------------------------------------------------------------
    
    
      Updated_Phe_Data[,metabolite_col_number] = log10( as.numeric( Updated_Phe_Data[,metabolite_col_number] ) )
      
      # # ----------------------------------------------------------------
    
    
      # # scaling CG values
      # # -----------------
    
      for ( c_loop in 361:ncol(Updated_Phe_Data) ) {
      
        Updated_Phe_Data[,c_loop] = scale( Updated_Phe_Data[,c_loop], center = TRUE, scale = TRUE )
      
      }
    
      # # ------------------------------------------------------------------------
    
    
      # independent variables in the regression
    
      Pheno_Header_List = c();                        Pheno_Header_List = names(Batch1_Phe)
    
      Independent_Variables_List = c();               Independent_Variables_List = Pheno_Header_List[96:131]  # CD8T to PC30_cp
    
      # constructing the right hand side of the regression equation with independent variables
      # right hand side of the regression equation contains CpG, age, gender, 6 cell types and 30 principal components associated with probes
    
      RHS_equation = '+Age+Sex'
    
      for (element_loop in Independent_Variables_List) {
      
        RHS_equation = paste(RHS_equation, '+', element_loop, sep="")
      
      }
    
    
    
      Results_bValues = c()
    
      Results_pValues = c()
    
      Results_seValues = c()
      
      Results_nValues = c()
    
      for ( CG_marker_loop in 1:length(Independent_Variables_List_CG) ) {  # 1:length(Independent_Variables_List_CG)
      
        NA_Indices = c();               NA_Indices = which( is.na(Updated_Phe_Data[,(360+CG_marker_loop)]) )
        
        # this 15% missing data filter is not actually used, as we already did sample call rate filtering
      
        if( length(NA_Indices) > round(0.15 * total_samples) ) {   # if more than 15% of the CG data is missing, we exclude it 
        
          Results_bValues = append( Results_bValues, 999 )
          
          Results_seValues = append( Results_seValues, 999 )
          
          Results_pValues = append( Results_pValues, 999 ) 
          
          Results_nValues = append( Results_nValues, 999 )
        
        } else {
      
          CG_marker_name = '';            CG_marker_name = Independent_Variables_List_CG[CG_marker_loop]
      
          # # set correct depVariable
          # # -----------------------
      
          lm_formula = '';                lm_formula = as.formula( as.character( paste(curr_metab_name, '~', CG_marker_name, RHS_equation, sep="") ) )    
                                                                                                                                                          
      
          lm_fit = lm(lm_formula, data = Updated_Phe_Data, na.action = na.exclude)
      
          Results_bValues = append( Results_bValues, coef( summary(lm_fit) )[2,1] )
          
          Results_seValues = append( Results_seValues, coef( summary(lm_fit) )[2,2] )
          
          Results_pValues = append( Results_pValues, coef( summary(lm_fit) )[2,4] ) 
          
          Results_nValues = append( Results_nValues, nrow(Updated_Phe_Data) ) 
      
        } # if-else loop
      
      } # CG_marker_loop
    
    
      LR_Results = data.frame( "CG_Marker" = Independent_Variables_List_CG, "bValues" = Results_bValues, "seValues" = Results_seValues, "pValues" = Results_pValues, "nValues" = Results_nValues )
    
      write.csv(LR_Results, file = results_filename, quote = FALSE, row.names = FALSE)
    
    
    } # job_number
  
    t2 = proc.time()
  
    print(paste('time taken for ', total_jobs, ' loops: ', sep=''))
  
    t2-t1
    
    rm(list = c('IA_epimigrant_beta','Updated_Phe_Data', 'Independent_Variables_Data'))
  
  
  } # batch tag loop
  
  