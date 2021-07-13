# Code description: 
# ------------------
# This code generates randomly sampled sets of Sentinel CpG markers and matching Background CpG markers. 
# For a user specified univariate association results data ( single phenotype vs. CpG levels ), this code performs permutation analysis.
# permutation analysis can be performed for any one of the three categoires of CpGs (a) cis, (b) long-range and (c) trans
# User can adjust (a) number of randomly sampled sets and (b) number of CpGs per sampled set. By default, 1000 sets are sampled, each with 1000 CpGs in it.

# Input files and user specified parameters
# ------------------------------------------
# (1) category tag: 'cis' or 'long-range' or 'trans'
# (2) phenotype name 
# (3) univariate results data file for the chosen phenotype
# (4) file containing the full list of Sentinel CpGs. This file should have a column specifying whether the CpG is 'cis' or 'long-range' or 'trans'
# (5) file containing the full list of cosmopolitan CpGs
# (6) file containing mean and SD values for each CpG marker


# Steps included in the code
# --------------------------
# (1) read the full Sentinel CpG list and extract separately the correct category of CpGs (i.e. 'cis' or 'long-range' or 'trans')
# (2) read list of cosmopolitan CpGs and generate Background CpGs set from Cosmopolitan set
# (3) generate randomly sampled sets of Sentinel and Baackground CpGs
# (4) perform permutation analysis
# (5) write permutation analysis results to output files. Save randomly sampled CpG sets to separate files

# Output files
# ------------
# (1) generated random samples of Sentinel and Background CpGs are saved as text files
# (2) permutation analysis results are saved as .csv files

# code last modified: 16/Nov/2020
# comments added later, without modifying the code 


  cat('\f')
  
  rm(list=ls())
  
  graphics.off()
  
  t0 = Sys.time()
  
  
  # # user specified data
  # # -------------------
  
  category_tag = 'cis'  # 'cis', 'long-range', 'trans'
  
  trait_name = 'BMI'
  
  
  # # Sentinel data
  # # -------------
  
  Sentinel_DF = read.table(file='/home/projects/12000713/Lakshmi/Perm_Testing_15Oct2020/st9_sentinel_snps_and_cpgs.txt', header = TRUE)
  
 
  # extracting the CpGs that belong to the user specified category 
  
  Category_Indices = c();               Category_Indices = which( Sentinel_DF$category == category_tag )
  
  Raw_Sentinel_CGs = c();               Raw_Sentinel_CGs = as.character( Sentinel_DF$cpg.sentinel[Category_Indices] )
  
  Sentinel_CGs = c();                   Sentinel_CGs = unique(Raw_Sentinel_CGs)
  
  sentinel_set_size = 0;                sentinel_set_size = min( 1000, length(Sentinel_CGs) )  # total sentinel CGs for long-range category is 499, hence we use, min(1000, 499) = 499  
  
  
  
  # # Cosmopolitan CGs List
  # # ----------------------
  
  load(file='/home/projects/12000713/Lakshmi/Perm_Testing_15Oct2020/Cosmopolitan_CGs.RData')
  
  NR_Cosmopolitan_CGs = c();          NR_Cosmopolitan_CGs = unique(Cosmopolitan_CGs_List)
  
  
  
  
  # # Trait univariate results data, for total list of CGs for which univariate association data is available
  # # -------------------------------------------------------------------------------------------------------
  
  load(file= paste('/home/projects/12000713/Lakshmi/Perm_Testing_15Oct2020/', trait_name, '_CG_B1_Overall_Results.RData', sep='') )  # Overall_Results
  
 
  Total_CGs_List = c();               Total_CGs_List = as.character( Overall_Results$CG_Marker )
   
  CGs_to_Filter = c();                CGs_to_Filter = setdiff(NR_Cosmopolitan_CGs, Total_CGs_List)          # CGs for which univariate results are not available
  
  Processed_Cosmo_CGs = c();          Processed_Cosmo_CGs = setdiff(NR_Cosmopolitan_CGs, CGs_to_Filter)
  
  Background_CGs = c();               Background_CGs = setdiff(Total_CGs_List, Processed_Cosmo_CGs)
  
  
  
  
  # # loading mean and SD values of CG markers in batch 1
  # # -----------------------------------------------------
  
  load(file='/home/projects/12000713/Lakshmi/Perm_Testing_15Oct2020/Batch1_CGs_Stats.RData') # CG_Stats_DF
  
  
  
  Temp_DF = data.frame('CG_Name' = Background_CGs)
  
  Background_CG_Stats_DF = merge(Temp_DF, CG_Stats_DF, by.x = 'CG_Name', all = FALSE, sort = FALSE)
  
 
  
  # Permutation analysis
  # --------------------
  
  Combined_Results = c()
  
  set.seed(1000000)    # for reproducibility of results, we need to use a same initial seed value 
  
 
  
  for(random_iteration in 1:1000) {
    
    
    Init_Sentinel_Set = c();                 Init_Sentinel_Set = sample(Sentinel_CGs, size = ( sentinel_set_size + round(0.1*sentinel_set_size) ), replace = FALSE )
    
    Init_Temp_DF = data.frame('CG_Name' = Init_Sentinel_Set)
    
    Init_Sent_CG_Stats_DF = merge(Init_Temp_DF, CG_Stats_DF, by.x = 'CG_Name', all = FALSE, sort = FALSE )
    
    
    Init_Sent_CG_Stats_DF$Mean_Lower_Cutoff = as.numeric(Init_Sent_CG_Stats_DF$Mean_Values) - ( (2/100) ) 
    
    Init_Sent_CG_Stats_DF$Mean_Upper_Cutoff = as.numeric(Init_Sent_CG_Stats_DF$Mean_Values) + ( (2/100) ) 
    
    
    Init_Sent_CG_Stats_DF$SD_Lower_Cutoff = as.numeric(Init_Sent_CG_Stats_DF$SD_Values) - ( 0.2 / 100 ) 
    
    Init_Sent_CG_Stats_DF$SD_Upper_Cutoff = as.numeric(Init_Sent_CG_Stats_DF$SD_Values) + ( 0.2 / 100 ) 
    
    
    Sampled_Sentinel_Set = c()
    
    # # generating the random background sample
    # # ---------------------------------------
    
    Temp_Background_CG_Stats = c();          Temp_Background_CG_Stats = Background_CG_Stats_DF   
    
    Sampled_Background_Set = c()
    
    successful_sample_count = 0;
    
    sentinel_entry_count = 0;
    
    
    while( (successful_sample_count < sentinel_set_size) & (sentinel_entry_count <= length(Init_Sentinel_Set)) )  {      # # | (sentinel_entry_count <= length(Init_Sentinel_Set)) 
      
      # print(sentinel_entry_count)
      
      sentinel_entry_count = sentinel_entry_count + 1
      
      current_sentinal_CG = '';               current_sentinal_CG =  as.character( Init_Sent_CG_Stats_DF$CG_Name[sentinel_entry_count] )   
      
      mean_lower_cutoff = 0;                 mean_lower_cutoff = Init_Sent_CG_Stats_DF$Mean_Lower_Cutoff[sentinel_entry_count]   
      
      mean_upper_cutoff = 0;                 mean_upper_cutoff = Init_Sent_CG_Stats_DF$Mean_Upper_Cutoff[sentinel_entry_count]    
      
      SD_lower_cutoff = 0;                   SD_lower_cutoff = Init_Sent_CG_Stats_DF$SD_Lower_Cutoff[sentinel_entry_count]    
      
      SD_upper_cutoff = 0;                   SD_upper_cutoff = Init_Sent_CG_Stats_DF$SD_Upper_Cutoff[sentinel_entry_count]    
      
      Compatible_Indices = c();
      
      Compatible_Indices = which( (Temp_Background_CG_Stats$Mean_Values >= mean_lower_cutoff) & 
                                  (Temp_Background_CG_Stats$Mean_Values <= mean_upper_cutoff) & 
                                  (Temp_Background_CG_Stats$SD_Values >= SD_lower_cutoff) & 
                                  (Temp_Background_CG_Stats$SD_Values <= SD_upper_cutoff) )
      
      
      
      if(length(Compatible_Indices) > 0 ) {
        
        Sampled_Sentinel_Set = append( Sampled_Sentinel_Set, current_sentinal_CG )
        
        chosen_index_value = 0;                chosen_index_value = sample(Compatible_Indices, 1, replace = FALSE)
        
        chosen_background_CG = '';             chosen_background_CG = as.character( Temp_Background_CG_Stats$CG_Name[chosen_index_value] )
      
        Sampled_Background_Set = append(Sampled_Background_Set, chosen_background_CG)
      
        # # removing the CG already sampled from background set
      
        Temp_Background_CG_Stats = Temp_Background_CG_Stats[-c(chosen_index_value),]
        
        successful_sample_count = successful_sample_count + 1
      
      } 
      
      
      
    } # while loop
    
    
    write.table(Sampled_Sentinel_Set, file = paste(category_tag,  '_Sentinel_CGs_iter', as.character(random_iteration), '.txt', sep=''), row.names = FALSE, col.names = FALSE, quote = FALSE )
    
    write.table(Sampled_Background_Set, file = paste(category_tag,  '_Background_CGs_iter', as.character(random_iteration), '.txt', sep=''), row.names = FALSE, col.names = FALSE, quote = FALSE )
    
    # # calculating the number of CGs associated with the trait, for a chosen p-value cutoff
    # # -------------------------------------------------------------------------------------
    
    pValues_List = c(0.05, 0.00005)
    
    Results_Data = c()
    
    Sentinel_Temp_DF = data.frame("CG_Marker" = Sampled_Sentinel_Set)
    
    Sentinel_EWAS_DF = merge(Sentinel_Temp_DF, Overall_Results, by.x = 'CG_Marker', all = FALSE, sort = FALSE )
    
    
    Background_Temp_DF = data.frame("CG_Marker" = Sampled_Background_Set)
    
    Background_EWAS_DF = merge(Background_Temp_DF, Overall_Results, by.x = 'CG_Marker', all = FALSE, sort = FALSE )
    
    
    for( pvalue_loop in 1:length(pValues_List) ){
      
      p_value = 99;                          p_value = pValues_List[pvalue_loop]
    
      sentinel_positive_hits = 0;            sentinel_positive_hits = length( which(Sentinel_EWAS_DF$pValues < p_value) )
      
      Results_Data = append(Results_Data, sentinel_positive_hits)
      
      background_positive_hits = 0;          background_positive_hits = length( which(Background_EWAS_DF$pValues < p_value) )
      
      Results_Data = append(Results_Data, background_positive_hits)
      
    }
    
    Results_Data = append(Results_Data, length(Sampled_Sentinel_Set))
    
    write.csv(Results_Data, file = paste( trait_name, '_', category_tag,  '_Permutation_Results_iter', as.character(random_iteration), '.csv', sep='' ), quote = FALSE, row.names = FALSE)
    
    Combined_Results = rbind(Combined_Results, Results_Data)
    
    
    
  } # random iteration loop
  
  t12 = Sys.time()
  
  colnames(Combined_Results) = c('Sent 0.05', 'Background 0.05', 'Sent 1E-5', 'Background 1E-5', 'set size')
  
  write.csv(Combined_Results, file = paste( trait_name, '_', category_tag,  '_Permutation_Results_Combined.csv', sep='' ), quote = FALSE, row.names = FALSE)
  
  print( 'time for the whole code' )
  
  print( t12 - t0 )

  
  
